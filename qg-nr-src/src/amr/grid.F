cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_point_contained:                                                     cc
cc                  Determines whether a given coordinate location            cc
cc                  occurs within a given grid.                               cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function grid_point_contained(gi,x,y,z)
      implicit none
      integer  gi
      real(kind=8) x, y, z
      include 'grid.inc'
      !include 'grid_methods.inc'
      !include 'glob.inc'
      logical     ltrace
      parameter ( ltrace = .false.)

      grid_point_contained = .false.

      if (  x .ge. gr_minx(gi)  .and.
     *      x .le. gr_maxx(gi)  .and.
     *      y .ge. gr_miny(gi)  .and.
     *      y .le. gr_maxy(gi)  .and.
     *      z .ge. gr_minz(gi)  .and.
     *      z .le. gr_maxz(gi)        ) then
         grid_point_contained = .true.
      end if

      if (ltrace) then
         write(*,*)'grid_point_contained: x/y/z:     ',x,y,z
         write(*,*)'grid_point_contained: minx/y/z:  ',
     *         gr_minx(gi),gr_miny(gi),gr_minz(gi)
         write(*,*)'grid_point_contained: maxx/y/z:  ',
     *         gr_maxx(gi),gr_maxy(gi),gr_maxz(gi)
         write(*,*)'grid_point_contained: ',grid_point_contained
      end if

      return
      end    ! END: grid_point_contained

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_point_is_refinedh:                                                   cc
cc                  Determines whether a given coordinate location            cc
cc                  occurs on a finer level with respect to the               cc
cc                  given grid.                                               cc
cc              NB: this version takes a grid resolution instead of           cc
cc                  grid number and is thereby able to be called by the       cc
cc                  hyper side of the code.                                   cc
cc        Returns 1 if this point occurs on a grid with finer resolution      cc
cc                0 if this point does not appear on a finer grid             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function grid_point_is_refinedh(x,y,z,myh)
      implicit none
      real(kind=8) x, y, z, myh
      !include 'grid_methods.inc'
      !include 'largesmall.inc'
      real(kind=8)        SMALLNUMBER
      parameter     (     SMALLNUMBER    =  1.0d-12  )
      integer     li, g2
      real(kind=8)grid_return_h
      external    grid_return_h
      integer     grid_return_sibling, level_return_start
      logical     grid_point_contained, grid_return_existence
      external    grid_point_contained, grid_return_existence
      logical     level_return_existence
      external    level_return_existence, level_return_start
      real(kind=8) hl
      integer     myid, proc_return_myid
      external          proc_return_myid
      logical     ltrace
      parameter ( ltrace = .false.)
      logical     ltrace2
      parameter ( ltrace2 = .false.)

      grid_point_is_refinedh = 0

      myid = proc_return_myid()
 98   format('[',I3,'] ',A,4G9.1)
 99   format('[',I3,'] ',A,2I5)
      if(ltrace2)write(*,77)'x/y/z/h', x,y,z,myh
      !if(ltrace)write(*,99)myid,' x/y/z/h', x,y,z,myh
 77   format(A,4G10.2)
 78   format(A,I3)
      !
      ! Only need to check next level
      ! since we are guarantee further levels of
      ! refinement must themselves be contained by
      ! the coarser levels:
      !
 
      ! First find first level finer than h:
      li = 1
  7   if (level_return_existence(li)) then
         g2 = level_return_start(li)
         hl = grid_return_h(g2)
         if (hl.gt. myh-SMALLNUMBER) then
            if(ltrace2)write(*,77)'  Coarse/Equal: hl,myh',hl,myh
            ! resolution is the same or coarser, look to next one:
            li = li + 1
            goto 7
         else
            ! this level is finer than h
            if(ltrace2)write(*,77)'  Finer         hl,myh',hl,myh
            goto 8
         end if
      else
         ! no level finer than h exists
         !grid_point_is_refinedh = 0
         if(ltrace2)write(*,77)'grid_point_is_refinedh:No finer level'
         !return
         goto 20
      end if

  8   continue
      if(ltrace2)write(*,78)'grid_point_is_refinedh:Check level: ',li

      if (level_return_existence(li)) then
         g2 = level_return_start(li)
         if(ltrace2)write(*,78)'grid_point_is_refinedh:Start at g2:',g2
  10     if (grid_return_existence(g2)) then
            if (grid_point_contained(g2,x,y,z)) then
               if(ltrace2)write(*,*)'grid_point_is_refinedh:Pt in g2',g2
               grid_point_is_refinedh = 1
               goto 20
            else
               if(ltrace2)write(*,*)'grid_point_is_refinedh:Pt notg2',g2
               g2 = grid_return_sibling(g2)
               goto 10
            end if
         end if
      end if
      
  20  continue

      if (ltrace2) then
         !write(*,*)'grid_point_is_refinedh:h:   '
         write(*,98)myid,'grid_point_is_refinedh:h:   ',myh
         write(*,98)myid,'grid_point_is_refinedh:x,y,z',x,y,z
         write(*,99)myid,'grid_point_is_refinedh',grid_point_is_refinedh
      else if (ltrace) then
         write(*,76)'x/y/z/h', x,y,z,myh,grid_point_is_refinedh
 76   format(A,4G10.2,I2)
      end if

      return
      end    ! END: grid_point_is_refinedh

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_point_is_refined:                                                    cc
cc                  Determines whether a given coordinate location            cc
cc                  occurs on a finer level with respect to the               cc
cc                  given grid.                                               cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function grid_point_is_refined(gi,x,y,z)
      implicit none
      integer  gi
      real(kind=8) x, y, z
      include 'grid_methods.inc'
      integer     li, g2
      logical     grid_point_contained
      external    grid_point_contained
      logical     ltrace
      parameter ( ltrace = .false.)

      grid_point_is_refined = .false.

      !
      ! Only need to check next level
      ! since we are guarantee further levels of
      ! refinement must themselves be contained by
      ! the coarser levels:
      !
      li = grid_return_level(gi) + 1
      if (level_return_existence(li)) then
         g2 = level_return_start(li)
  10     if (grid_return_existence(g2)) then
            if (grid_point_contained(g2,x,y,z)) then
               if(ltrace)write(*,*)'grid_point_is_refined:Pt in  g2:',g2
               grid_point_is_refined = .true.
               goto 20
            else
               if(ltrace)write(*,*)'grid_point_is_refined:Pt not g2:',g2
               g2 = grid_return_sibling(g2)
               goto 10
            end if
         end if
      end if
      
  20  continue

      if (ltrace) then
         write(*,*)'grid_point_is_refined: gi:   ',gi
         write(*,*)'grid_point_is_refined: x,y,z:',x,y,z
         write(*,*)'grid_point_is_refined:       ',grid_point_is_refined
      end if

      return
      end    ! END: grid_point_is_refined

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  load_pointers:                                                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine load_pointers(grid)
      implicit none
      integer  grid
      include 'grid.inc'
      integer  i

      do i = 1, num_gfuncs
         gfunc_tmpptr(i) = gfunc_pointer(i,grid)
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_set_mask:                                                            cc
cc                  Called such that this process owns g1 and sets            cc
cc                  the chr() mask function only on g1 where g1's             cc
cc                  boundary lies in interior of g2.                          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_set_mask(g1,g2)
      implicit none
      integer  g1, g2
      include 'grid.inc'
      include 'grid_methods.inc'
      include 'glob.inc'
      integer     myid
      integer     proc_return_myid
      external    proc_return_myid
      integer  length, min_i, max_i, min_io, max_io,
     *                 min_j, max_j, min_jo, max_jo,
     *                 min_k, max_k, min_ko, max_ko,
     *                 nx, ny, nz,   nxo, nyo, nzo,
     *                 smallx, largex,
     *                 smally, largey,
     *                 smallz, largez, minspan
      logical     double_equal
      external    double_equal
      character(3) iterch1, iterch2
      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )

      myid = proc_return_myid()
 99   format('[',I3,'] ',A,2I5)
      if (ltrace) then
         write(*,99)myid, 'grid_set_mask: g1, g2 = ', g1, g2
         if (ltrace2) call grid_dump_info(g1)
      end if

      if (.not.grid_is_local(g1)) then
         write(*,99)myid,'grid_set_mask:Problem g1 is not locally owned'
         call my_exit('grid_set_mask: Problem, g1 is not locally owned')
      end if

      call grid_find_intersection(g1, g2,length,
     *            min_i, max_i, min_j, max_j, min_k, max_k)

      if (length.eq.0) then
         if(ltrace)write(*,99)myid,'   No intersection. Returning.'
         return
      end if

      !
      ! Check that grids overlap sufficiently so that boundary
      ! data of one is in the interior of other:
      !
      minspan = min(max_i-min_i+1, max_j-min_j+1, max_k-min_k+1)
      if (minspan .le. bound_width) then
         ! This is important, for example, when two refined regions
         ! in orbit about each other, overlap just a couple points
         ! if we were to overwrite the AMR_BDY values, then the
         ! boundary points would not get set by syncing or by interpolation from parents
         if(ltrace)write(*,99)myid,'grid_set_mask:Insufficientintersect'
         return
      end if

      call grid_find_intersection(g2, g1, length,
     *            min_io,max_io,min_jo,max_jo,min_ko,max_ko)

      if(ltrace2)write(*,*)myid,'grid_set_mask: intersections computed'

      nx  = gr_nx(g1)
      ny  = gr_ny(g1)
      nz  = gr_nz(g1)
      !
      nxo = gr_nx(g2)
      nyo = gr_ny(g2)
      nzo = gr_nz(g2)

      if (ltrace) then
       write(*,*)'grid_set_mask: nx,ny,nz    = ',nx,ny,nz
       write(*,*)'grid_set_mask: gr_min/maxx = ',gr_minx(g1),gr_maxx(g1)
       write(*,*)'grid_set_mask: gr_min/maxy = ',gr_miny(g1),gr_maxy(g1)
       write(*,*)'grid_set_mask: gr_min/maxz = ',gr_minz(g1),gr_maxz(g1)
       write(*,*)'grid_set_mask: min/max x 0 = ',minx0,maxx0
       write(*,*)'grid_set_mask: min/max y 0 = ',miny0,maxy0
       write(*,*)'grid_set_mask: min/max z 0 = ',minz0,maxz0
       !write(*,*)'grid_set_mask: ',gr_maxx(g1).eq.maxx0
       !write(*,*)'grid_set_mask: ',gr_maxy(g1).eq.maxy0
       !write(*,*)'grid_set_mask: ',gr_maxz(g1).eq.maxz0
       !write(*,*)'grid_set_mask: ',gr_maxx(g1)-maxx0
       !write(*,*)'grid_set_mask: ',gr_maxy(g1)-maxy0
       !write(*,*)'grid_set_mask: ',gr_maxz(g1)-maxz0
       !write(*,*)'grid_set_mask: small/largex= ',smallx,largex
       !write(*,*)'grid_set_mask: small/largey= ',smally,largey
       !write(*,*)'grid_set_mask: small/largez= ',smallz,largez
      end if

      call load_pointers(g1)

      if (ltrace2) write(*,*) 'grid_set_mask: Call setmsk_int_bnd'
      if (ltrace2) then
         write(*,*) '      min_i,   max_i   = ', min_i,   max_i
         write(*,*) '      min_io,  max_io  = ', min_io,  max_io
      end if

      call setmsk_int_bndc(   q(gr_chr),
     *            min_i, max_i, min_j, max_j, min_k, max_k,
     *            min_io,max_io,min_jo,max_jo,min_ko,max_ko,
     *            nx,ny,nz, nxo,nyo,nzo, ghostwidth,bound_width)
!     call setmsk_int_bndB(   q(gr_chr),
!    *            min_i, max_i, min_j, max_j, min_k, max_k,
!    *            min_io,max_io,min_jo,max_jo,min_ko,max_ko,
!    *            smallx, largex, 
!    *            smally, largey, 
!    *            smallz, largez, 
!    *            nx,ny,nz, nxo,nyo,nzo, ghostwidth,bound_width)

      if (ltrace2) then
         call grid_dump_info(g1)
         write(*,99)myid,'grid_set_mask: Done.',g1,g2
                call int2str(g1,iterch1)
                call int2str(g2,iterch2)
               call field_out3d(q(gr_chr),
     *             gr_t(g1),'mask_chr'//iterch1//iterch2,
     *             gr_minx(g1),gr_maxx(g1),
     *             gr_miny(g1),gr_maxy(g1),
     *             gr_minz(g1),gr_maxz(g1),
     *             gr_nx(g1),gr_ny(g1),gr_nz(g1),myid)
      end if

      return
      end       ! END: grid_set_mask

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_set_amrmask:                                                         cc
cc                  Go back and ensure AMR boundaries are well marked.        cc
cc                  (1) Fix grid edges which have been set as DECO before.    cc
cc                  (2) Fix lines of intersection of DECO on faces:           cc
cc                  (3) Ensure AMR boundary is backed by interior points.     cc
cc                  These tasks must happen after all grids are looked at to  cc
cc                  set the DECO boundary.                                    cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_set_amrmask(gi)
      implicit none
      integer  gi
      include 'grid.inc'
      include 'grid_methods.inc'
      include 'glob.inc'
      include 'chr.inc'
      integer     myid
      integer     proc_return_myid
      external    proc_return_myid
      integer     i,j,k, index, indexjp1, indexkp1, ll,
     *                 nx, ny, nz, ii,jj,kk, indexii,indexjj,
     *            previ_pt, prevj_pt, cur_pt, tmpi
      integer     AMR_BDY, DECOBDY, INTERIOR, neighbors
      logical     nearamrbdy
      logical     ltrace
      character(3) iterch
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )

      !ltrace  = .false.
      !ltrace2 = .false.
      if (.false.) then
         write(*,99)myid,'grid_set_amrmask:Bypassing and returning'
         return
      end if 

      AMR_BDY  = NINT(CHR_AMR_BDY)
      DECOBDY  = NINT(CHR_DECO_BDY)

      myid = proc_return_myid()
 99   format('[',I3,'] ',A,2I5)
      if (ltrace) then
         write(*,99)myid, 'grid_set_amrmask: gi=',gi
         if (ltrace2) call grid_dump_info(gi)
      end if

      if (.not.grid_is_local(gi)) then
         write(*,99)myid,'grid_set_amrmask:Problem gi not locally owned'
         call my_exit('grid_set_amrmask: Problem, gi not locally owned')
      end if

      nx  = gr_nx(gi)
      ny  = gr_ny(gi)
      nz  = gr_nz(gi)
      !if (nx.eq.53.and.ny.eq.53.and.nz.eq.53) then
         !ltrace  = .true.
         !ltrace2 = .true.
      !end if
      !
      call load_pointers(gi)

      ! Work on copy of chr() in order to "grow" amr regions a bit:
      call mat_copy3d(q(gr_chr), q(gr_tmp), nx,ny,nz)
      !
         if (ltrace2) write(*,99)myid,'grid_setamrmask: The alternative'
         ! Alternative way:
         do k = 1, nz, nz-1
            !Z-faces
            do j = 1, ny
            do i = 1, nx
            !do j = 1+bound_width-1, ny-bound_width+1
            !do i = 1+bound_width-1, nx-bound_width+1
               index    = (k-1)*nx*ny+(j-1)*nx+(i-1)
               if (NINT(q(gr_chr+index)).eq.DECOBDY) then
                 nearamrbdy = .false.
                 ! we do not need to consider the corners:
                 do jj=max(j-1,2), min(j+1,ny-1)
                 indexii  = (k-1)*nx*ny+(jj-1)*nx+(i-1)
                 if(NINT(q(gr_chr+indexii)).eq.AMR_BDY)nearamrbdy=.true.
                 end do
                 do ii=max(i-1,2), min(i+1,nx-1)
                 indexii  = (k-1)*nx*ny+(j-1)*nx+(ii-1)
                 if(NINT(q(gr_chr+indexii)).eq.AMR_BDY)nearamrbdy=.true.
                 end do
                 if (nearamrbdy) then
                    ! Change points directly behind this point to interior:
                    if(ltrace)write(*,*)'grid_set_amrmask: Pt:', i, j, k
                    do kk=max(k-bound_width+1,1),min(k+bound_width-1,nz)
                       ! Make an expanding triangle behind point interior:
                       ! assymmetry is so that each point is interior in
                       ! *only* one of the overlapping grids:
                       tmpi = abs(kk-k)
                       !if(ltrace)write(*,*)'grid_set_amrmask:tm',tmpi,kk
                       do jj= max(j-tmpi,   bound_width+2-tmpi),
     .                        min(j+tmpi-1, ny-bound_width+tmpi-1)
                       do ii= max(i-tmpi,   bound_width+2-tmpi),
     .                        min(i+tmpi-1, nx-bound_width+tmpi-1)
                       !do ii=max(i-tmpi,2),min(i+tmpi-1,nx-1)
                       !do jj=max(j-tmpi,2),min(j+tmpi-1,ny-1)
                          indexii  = (kk-1)*nx*ny+(jj-1)*nx+(ii-1)
                          if(ltrace)write(*,*)'amrmask:',ii,jj,kk
                          indexjj  = (k-1)*nx*ny+(jj-1)*nx+(ii-1)
                          if ( (ii.eq.i.and.jj.eq.j) .or.
     .                         (NINT(q(gr_chr+indexjj)).eq.DECOBDY))then
                             ! only expand triangle directly behind
                             ! point and towards the region that is
                             ! not an AMR boundary:
                             q(gr_tmp+indexii) = CHR_INTERIOR
                          end if
                       end do
                       end do
                    end do
                    ! Certain points on the boundary are missed by other code,
                    ! So this ensures it is set correctly:
                    q(gr_tmp+index) = CHR_AMR_BDY
                 end if
               end if
            end do
            end do
         end do
         if (ltrace) write(*,99)myid,'grid_setamrmask: Ybounds:'
         do j = 1, ny, ny-1
            !Y-faces
            do k = 1, nz
            do i = 1, nx
            !do k = 1+bound_width-1, nz-bound_width+1
            !do i = 1+bound_width-1, nx-bound_width+1
               index    = (k-1)*nx*ny+(j-1)*nx+(i-1)
               if (NINT(q(gr_chr+index)).eq.DECOBDY) then
                 nearamrbdy = .false.
                 do kk=max(k-1,2), min(k+1,nz-1)
                 indexii  = (kk-1)*nx*ny+(j-1)*nx+(i-1)
                 if(NINT(q(gr_chr+indexii)).eq.AMR_BDY)nearamrbdy=.true.
                 end do
                 do ii=max(i-1,2), min(i+1,nx-1)
                 indexii  = (k-1)*nx*ny+(j-1)*nx+(ii-1)
                 if(NINT(q(gr_chr+indexii)).eq.AMR_BDY)nearamrbdy=.true.
                 end do
                 if (nearamrbdy) then
                    ! Change points directly behind this point to interior:
                    if(ltrace)write(*,*)'grid_set_amrmask:Ybnd Pt',i,j,k
                    do jj=max(j-bound_width+1,1),min(j+bound_width-1,ny)
                       ! Make an expanding triangle behind point interior:
                       ! a symmetry is so that each point is interior in
                       ! *only* one of the overlapping grids:
                       tmpi = abs(jj-j)
                       !if(ltrace)write(*,*)'grid_set_amrmask:tm',tmpi,jj
                       !do ii=max(i-tmpi,2),min(i+tmpi-1,nx-1)
                       !do kk=max(k-tmpi,2),min(k+tmpi-1,nz-1)
                       do kk= max(k-tmpi,   bound_width+2-tmpi),
     .                        min(k+tmpi-1, nz-bound_width+tmpi-1)
                       do ii= max(i-tmpi,   bound_width+2-tmpi),
     .                        min(i+tmpi-1, nx-bound_width+tmpi-1)
                          indexii  = (kk-1)*nx*ny+(jj-1)*nx+(ii-1)
                          if(ltrace)write(*,*)'amrmask:',ii,jj,kk
                          indexjj  = (kk-1)*nx*ny+(j-1)*nx+(ii-1)
                          if ( (ii.eq.i.and.kk.eq.k) .or.
     .                         (NINT(q(gr_chr+indexjj)).eq.DECOBDY))then
                             ! only expand triangle directly behind
                             ! point and towards the region that is
                             ! not an AMR boundary:
                             q(gr_tmp+indexii) = CHR_INTERIOR
                          end if
                       end do
                       end do
                    end do
                    ! Certain points on the boundary are missed by other code,
                    ! So this ensures it is set correctly:
                    q(gr_tmp+index) = CHR_AMR_BDY
                 end if
               end if
            end do
            end do
         end do
      if (.true.) then
         do i = 1, nx, nx-1
            !X-faces
            do k = 1, nz
            do j = 1, ny
            !do k = 1+bound_width-1, nz-bound_width+1
            !do j = 1+bound_width-1, ny-bound_width+1
               index    = (k-1)*nx*ny+(j-1)*nx+(i-1)
               if (NINT(q(gr_chr+index)).eq.DECOBDY) then
                 nearamrbdy = .false.
                 do kk=max(k-1,2), min(k+1,nz-1)
                 indexii  = (kk-1)*nx*ny+(j-1)*nx+(i-1)
                 if(NINT(q(gr_chr+indexii)).eq.AMR_BDY)nearamrbdy=.true.
                 end do
                 do jj=max(j-1,2), min(j+1,ny-1)
                 indexii  = (k-1)*nx*ny+(jj-1)*nx+(i-1)
                 if(NINT(q(gr_chr+indexii)).eq.AMR_BDY)nearamrbdy=.true.
                 end do
                 if (nearamrbdy) then
                    ! Change points directly behind this point to interior:
                    if(ltrace)write(*,*)'grid_set_amrmask: Pt:', i, j, k
                    do ii=max(i-bound_width+1,1),min(i+bound_width-1,nx)
                       ! Make an expanding triangle behind point interior:
                       ! a symmetry is so that each point is interior in
                       ! *only* one of the overlapping grids:
                       tmpi = abs(ii-i)
                       if(ltrace)write(*,*)'grid_set_amrmask:tm',tmpi,ii
                       do kk= max(k-tmpi,   bound_width+2-tmpi),
     .                        min(k+tmpi-1, nz-bound_width+tmpi-1)
                       do jj= max(j-tmpi,   bound_width+2-tmpi),
     .                        min(j+tmpi-1, ny-bound_width+tmpi-1)
                       !do jj=max(j-tmpi,2),min(j+tmpi-1,ny-1)
                       !do kk=max(k-tmpi,2),min(k+tmpi-1,nz-1)
                          indexii  = (kk-1)*nx*ny+(jj-1)*nx+(ii-1)
                          if(ltrace)write(*,*)'amrmask:',ii,jj,kk
                          indexjj  = (kk-1)*nx*ny+(jj-1)*nx+(i-1)
                          if ( (jj.eq.j.and.kk.eq.k) .or.
     .                         (NINT(q(gr_chr+indexjj)).eq.DECOBDY))then
                             ! only expand triangle directly behind
                             ! point and towards the region that is
                             ! not an AMR boundary:
                             q(gr_tmp+indexii) = CHR_INTERIOR
                          end if
                       end do
                       end do
                    end do
                    ! Certain points on the boundary are missed by other code,
                    ! So this ensures it is set correctly:
                    q(gr_tmp+index) = CHR_AMR_BDY
                 end if
               end if
            end do
            end do
         end do
      end if

      ! the above algorithm can produce some strange looking
      ! "interference" pattern near the edges where two
      ! boundaries come together...this is an attempt to
      ! clean it up a bit...in particular small regions
      ! of interior points are created within a sea
      ! of domain decomposition points that can cause
      ! some nonsmoothness.
      !   And so here, look for any Interior points
      ! that were DD points but now have at least
      ! 4 dd neighbors:
      INTERIOR = NINT(CHR_INTERIOR)
      do ll = 1, bound_width
         do k = 2, nz-1
         do j = 2, ny-1
         do i = 2, nx-1
            indexii  = (k-1)*nx*ny+(j-1)*nx+(i-1)
            if (NINT(q(gr_tmp+indexii)).eq.INTERIOR .and.
     *          NINT(q(gr_chr+indexii)).eq.DECOBDY) then
               ! This is a point that has been
               ! "switched" in the above code:
               !
               ! Cound how many DD neighbors:
               neighbors = 0
               indexjj  = (k-1)*nx*ny+(j-1)*nx+(i-2)
               if (NINT(q(gr_tmp+indexjj)).eq.DECOBDY)
     .                           neighbors=neighbors+1
               indexjj  = (k-1)*nx*ny+(j-2)*nx+(i-1)
               if (NINT(q(gr_tmp+indexjj)).eq.DECOBDY)
     .                           neighbors=neighbors+1
               indexjj  = (k-2)*nx*ny+(j-1)*nx+(i-1)
               if (NINT(q(gr_tmp+indexjj)).eq.DECOBDY)
     .                           neighbors=neighbors+1
               indexjj  = (k-1)*nx*ny+(j-1)*nx+(i  )
               if (NINT(q(gr_tmp+indexjj)).eq.DECOBDY)
     .                           neighbors=neighbors+1
               indexjj  = (k-1)*nx*ny+(j  )*nx+(i-1)
               if (NINT(q(gr_tmp+indexjj)).eq.DECOBDY)
     .                           neighbors=neighbors+1
               indexjj  = (k  )*nx*ny+(j-1)*nx+(i-1)
               if (NINT(q(gr_tmp+indexjj)).eq.DECOBDY)
     .                           neighbors=neighbors+1
               !
               if (neighbors.ge.4) then
                  q(gr_tmp+indexii) = CHR_DECO_BDY
               end if
            end if
         end do
         end do
         end do
      end do

      if (ltrace) then
                call int2str(gi,iterch)
               call field_out3d(q(gr_chr),
     *             gr_t(gi),'amrmA_chr'//iterch,
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             gr_nx(gi),gr_ny(gi),gr_nz(gi),myid)
      end if

      ! Copy new chr() array back into place:
      call mat_copy3d(q(gr_tmp), q(gr_chr), nx,ny,nz)
      call load_scal3d(q(gr_tmp), 0.d0, nx,ny,nz )

      if (ltrace) then
               call field_out3d(q(gr_chr),
     *             gr_t(gi),'amrmB_chr'//iterch,
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             gr_nx(gi),gr_ny(gi),gr_nz(gi),myid)
      end if

      return
      end       ! END: grid_set_amrmask


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_dump_info:                                                           cc
cc                                                                            cc
cc      Output grid info.                                                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_dump_info(gi)
      implicit none
      integer  gi
      include 'grid.inc'
      include 'mpif.h'
      include 'mpi_stuff.inc'
      integer  i
      real(kind=8) myl2norm3d, mynorm
      external     myl2norm3d

 96   format('[',I3,']   Problem: ',I4,A,3G21.14)
 97   format('[',I3,']     ',I4,A,3F21.14)
 98   format('[',I3,'] ',A,3F15.8)
 99   format('[',I3,'] ',A,3I5)

      write(*,99) myid, ' ----------grid_dump_info:'
      write(*,99) myid, ' GRID:       ',gi
      write(*,99) myid, ' nx,ny,nz:   ',gr_nx(gi),gr_ny(gi),gr_nz(gi)
      write(*,99) myid, ' owner:      ',gr_own(gi)
      write(*,99) myid, ' level:      ',gr_level(gi)
      write(*,98) myid, ' h:          ',gr_h(gi)
      write(*,98) myid, ' t:          ',gr_t(gi)
      write(*,99) myid, ' gr_sibling: ',gr_sibling(gi)
      write(*,98) myid, ' minx/y/z:   ',gr_minx(gi),gr_miny(gi),
     *                                  gr_minz(gi)
      write(*,98) myid, ' maxx/y/z:   ',gr_maxx(gi),gr_maxy(gi),
     *                                  gr_maxz(gi)
      !
      ! If grid is local, we have more info to print:
      !
      if (gr_own(gi).eq.myid) then
          write(*,*) ' x(1),x(nx):',
     *                         q(gr_x(gi)),q(gr_x(gi)+gr_nx(gi)-1)
          write(*,*) ' y(1),y(ny):',
     *                         q(gr_y(gi)),q(gr_y(gi)+gr_ny(gi)-1)
          write(*,*) ' z(1),z(nz):',
     *                         q(gr_z(gi)),q(gr_z(gi)+gr_nz(gi)-1)
          write(*,*) ' gr_x,gr_y,gr_z:',gr_x(gi),gr_y(gi),gr_z(gi)
          do i = 1, num_gfuncs
            mynorm =
     * myl2norm3d(q(gfunc_pointer(i,gi)),gr_nx(gi),gr_ny(gi),gr_nz(gi))
            !if (mynorm.gt.1d8.or.i.eq.1)
            call field_dump_infob(q(gfunc_pointer(i,gi)),
     *           gr_nx(gi),gr_ny(gi),gr_nz(gi),gfunc_name(i)(1:10))
            !write(*,97) myid, i,gfunc_name(i)(1:7), mynorm
            if (mynorm.gt.1d8)
     *         write(*,96) myid, i,gfunc_name(i)(1:7), mynorm
!           call field_out3d(q(gfunc_pointer(i,gi)),gr_t(gi),
!    *           'dump',
!    *             gr_minx(gi),gr_maxx(gi),
!    *             gr_miny(gi),gr_maxy(gi),
!    *             gr_minz(gi),gr_maxz(gi),
!    *             gr_nx(gi),gr_ny(gi),gr_nz(gi),myid)
          end do
c            i = 61
c            write(*,97) myid, i,gfunc_name(i)(1:7),
c    * myl2norm3d(q(gfunc_pointer(i,gi)),gr_nx(gi),gr_ny(gi),gr_nz(gi))
      end if

      write(*,99) myid, ' ----------grid_dump_info:'

      return
      end       ! END: grid_dump_info

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_get_dims:                                                            cc
cc                                                                            cc
cc      Return the dimensions of the grid's domain                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_get_dims(grid,nx,ny,nz)
      implicit none
      include 'grid.inc'
      integer  grid, nx,ny,nz
      logical     ltrace
      parameter ( ltrace = .false. )

      if(ltrace)write(*,*)'grid_get_dims: Enter: ',grid
      nx = gr_nx(grid)
      ny = gr_ny(grid)
      nz = gr_nz(grid)
      if(ltrace)write(*,*)'grid_get_dims: Done.: ',grid

      return
      end       ! END: grid_get_dims

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_mk_pending:                                                          cc
cc                                                                            cc
cc      Set grid to have status Pending.                                      cc
cc      Upon entering the refinement process for a given level,               cc
cc      all pre-existing grids are marked DEAD. They cannot be deleted        cc
cc      yet because we need their data to initialize new grids. Plus,         cc
cc      some of these grids might get resuse in which case we can just        cc
cc      mark them alive.                                                      cc
cc      However, we need an intermediate step in which grids are marked       cc
cc      as Pending. If the refinement process completes without running       cc
cc      out of memory, then all pending grids will be marked as alive.        cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_mk_pending(grid)
      implicit none
      include 'grid.inc'
      integer  grid

      gr_alive(grid) = PENDING

      return
      end       ! END: grid_mk_pending

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_mk_dead:                                                             cc
cc                                                                            cc
cc      Set grid to be dead                                                   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_mk_dead(grid)
      implicit none
      include 'grid.inc'
      integer  grid

      gr_alive(grid) = DEAD

      return
      end       ! END: grid_mk_dead

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_mk_alive:                                                            cc
cc                                                                            cc
cc      Set grid to be alive:                                                 cc
cc                                                                            cc
cc                     gr_alive(grid) = ALIVE                                 cc
cc                                                                            cc
cc      Also re-initializes the mask function chr() because                   cc
cc      its overlapping boundaries will have changed and need                 cc
cc      to be determined anew.                                                cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_mk_alive(grid)
      implicit none
      include 'grid.inc'
      include 'grid_methods.inc'
      integer  grid

      logical     ltrace
      parameter ( ltrace = .false. )

      if (ltrace) write(*,*) 'grid_mk_alive: Enter, grid=',grid

      gr_alive(grid) = ALIVE

      call load_pointers(grid)

      if (ltrace) write(*,*) 'grid_mk_alive: Before grid_init_chr'

      if (grid_is_local(grid)) then
         if (ltrace) write(*,*) 'grid_mk_alive: Initting chr array'
         call grid_init_chr(   q(gr_chr), 
     *                      gr_minx(grid), gr_miny(grid),
     *                      gr_minz(grid), gr_h(grid),
     *                      gr_nx(grid),   gr_ny(grid),
     *                      gr_nz(grid)                      )
      end if

      if (ltrace) write(*,*) 'grid_mk_alive: Done.'

      return
      end       ! END: grid_mk_alive

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_set_parent:                                                          cc
cc                  Set the parent field of a grid to a specified grid.       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_set_parent( gi, newparent )
      implicit    none
      integer     gi, newparent
      include    'grid.inc'
      
      gr_parent( gi ) = newparent

      return
      end       ! END: grid_set_parent

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_set_child:                                                           cc
cc                  Set the child  field of a grid to a specified grid.       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_set_child( gi, newchild )
      implicit    none
      integer     gi, newchild
      include    'grid.inc'
      
      gr_child( gi ) = newchild

      return
      end       ! END: grid_set_child

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_set_sibling:                                                         cc
cc                  Set the sibling ield of a grid to a specified grid.       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_set_sibling( gi, newsibling )
      implicit    none
      integer     gi, newsibling
      include    'grid.inc'
      
      gr_sibling( gi ) = newsibling

      return
      end       ! END: grid_set_sibling

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_contained_within:                                                    cc
cc                  Is grid "gi" contained within grid "gp"?                  cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function grid_contained_within( gi, gp )
      implicit    none
      integer     gi, gp
      include    'grid.inc'
      logical     ltrace
      parameter ( ltrace = .false. )
      
      if (ltrace) then
          write(*,*) 'Possible Parent:'
          write(*,*) ' x: ',gr_minx(gp),gr_maxx(gp)
          write(*,*) ' y: ',gr_miny(gp),gr_maxy(gp)
          write(*,*) ' z: ',gr_minz(gp),gr_maxz(gp)
          write(*,*) 'Child:'
          write(*,*) ' x: ',gr_minx(gi),gr_maxx(gi)
          write(*,*) ' y: ',gr_miny(gi),gr_maxy(gi)
          write(*,*) ' z: ',gr_minz(gi),gr_maxz(gi)
      end if

      grid_contained_within =   gr_minx(gi) .ge. gr_minx(gp)
     *                    .and. gr_miny(gi) .ge. gr_miny(gp)
     *                    .and. gr_minz(gi) .ge. gr_minz(gp)
     *                    .and. gr_maxx(gi) .le. gr_maxx(gp)
     *                    .and. gr_maxy(gi) .le. gr_maxy(gp)
     *                    .and. gr_maxz(gi) .le. gr_maxz(gp)

      if (ltrace) then
         write(*,*) 'grid_contained_within: gi = ',gi
         write(*,*) 'grid_contained_within: gp = ',gp
         write(*,*) 'grid_contained_within:    =',grid_contained_within
      end if

      return
      end       ! END: grid_contained_within

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_return_work:                                                         cc
cc                  Returns work value associated to the input grid.          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8) function grid_return_work(gridnumber)
      implicit    none
      include    'grid.inc'
      integer     gridnumber

      grid_return_work = gr_work(gridnumber)

      return
      end       ! END: grid_return_work

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_intersect:                                                           cc
cc           Logical: true  if grids have nonzero intersection                cc
cc                    false if grids have no      interesction                cc
cc                                                                            cc
cc           Grids need not be on the same level.                             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function grid_intersect(g1, g2)
      implicit    none
      integer     g1, g2
      include    'grid.inc'
      include    'largesmall.inc'

      logical     ltrace
      parameter ( ltrace = .false.)

      if (ltrace) write(*,*) 'grid_intersect: g1,g2=',g1,g2

      grid_intersect = .true.

      if ( gr_minx(g1) .ge. gr_minx(g2)-SMALLNUMBER ) then
         if ( gr_minx(g1) .gt. gr_maxx(g2)+SMALLNUMBER ) then
            grid_intersect = .false.
            goto 100
         end if
      else
         if ( gr_maxx(g1) .lt. gr_minx(g2)-SMALLNUMBER ) then
            grid_intersect = .false.
            goto 100
         end if
      end if
      !
      if ( gr_miny(g1) .ge. gr_miny(g2)-SMALLNUMBER ) then
         if ( gr_miny(g1) .gt. gr_maxy(g2)+SMALLNUMBER) then
            grid_intersect = .false.
            goto 100
         end if
      else
         if ( gr_maxy(g1) .lt. gr_miny(g2)-SMALLNUMBER ) then
            grid_intersect = .false.
            goto 100
         end if
      end if
      !
      if ( gr_minz(g1) .ge. gr_minz(g2)-SMALLNUMBER ) then
         if ( gr_minz(g1) .gt. gr_maxz(g2)+SMALLNUMBER) then
            grid_intersect = .false.
            goto 100
         end if
      else
         if ( gr_maxz(g1) .lt. gr_minz(g2)-SMALLNUMBER ) then
            grid_intersect = .false.
            goto 100
         end if
      end if

 100  continue

      if (ltrace) then
          write(*,*) 'grid_intersect = ',grid_intersect
      end if

      return
      end       ! END: grid_intersect

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_find_intersection:                                                   cc
cc                  Returns the intersection between grid1 and grid2          cc
cc                  relative to grid1.                                        cc
cc                  Assumes grid on same level.                               cc
cc            Returns: length of vector needed to store intersection.         cc
cc                     the start (si,sj,sk) and finish (fi,fj,fk) of          cc
cc                        the intersection region within grid1.               cc
cc            NB: this one doesn't use anchors and therefore can deal         cc
cc                with grids from different parents.                          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_find_intersection(g1, g2, length, si,fi,
     *                  sj,fj, sk,fk)
      implicit    none
      integer     g1, g2
      integer     si,fi, sj,fj, sk,fk, length
      include    'grid.inc'
      include    'glob.inc'

      logical     ltrace
      parameter ( ltrace = .false.)

      if (ltrace) write(*,*) 'grid_find_intersection: g1,g2=',g1,g2

      length = 0
      si     = 0
      fi     = 0
      sj     = 0
      fj     = 0
      sk     = 0
      fk     = 0

      ! 
      ! Find bounding box of intersection (if any):
      !
      if ( gr_minx(g1) .ge. gr_minx(g2) ) then
         if ( gr_minx(g1) .gt. gr_maxx(g2) ) goto 100! no intersection
         si = 1
      else
         if ( gr_maxx(g1) .lt. gr_minx(g2) ) goto 100! no intersection
         si = 1 + NINT( (gr_minx(g2) - gr_minx(g1)) / gr_h(g1) )
      end if
      if ( gr_maxx(g1) .le. gr_maxx(g2) ) then
         fi = gr_nx(g1)
      else
         fi = 1 + NINT( (gr_maxx(g2) - gr_minx(g1)) / gr_h(g1) )
      end if
      !
      if ( gr_miny(g1) .ge. gr_miny(g2) ) then
         if ( gr_miny(g1) .gt. gr_maxy(g2) ) goto 100! no intersection
         sj = 1
      else
         if ( gr_maxy(g1) .lt. gr_miny(g2) ) goto 100! no intersection
         sj = 1 + NINT( (gr_miny(g2) - gr_miny(g1)) / gr_h(g1) )
      end if
      if ( gr_maxy(g1) .le. gr_maxy(g2) ) then
         fj = gr_ny(g1)
      else
         fj = 1 + NINT( (gr_maxy(g2) - gr_miny(g1)) / gr_h(g1) )
      end if
      !
      if ( gr_minz(g1) .ge. gr_minz(g2) ) then
         if ( gr_minz(g1) .gt. gr_maxz(g2) ) goto 100! no intersection
         sk = 1
      else
         if ( gr_maxz(g1) .lt. gr_minz(g2) ) goto 100! no intersection
         sk = 1 + NINT( (gr_minz(g2) - gr_minz(g1)) / gr_h(g1) )
      end if
      if ( gr_maxz(g1) .le. gr_maxz(g2) ) then
         fk = gr_nz(g1)
      else
         fk = 1 + NINT( (gr_maxz(g2) - gr_minz(g1)) / gr_h(g1) )
      end if

      length = (fk-sk+1)*(fj-sj+1)*(fi-si+1)

 100  continue

      if (ltrace) then
          write(*,*) '        minx1,maxx1 = ',gr_minx(g1),gr_maxx(g1)
          write(*,*) '        miny1,maxy1 = ',gr_miny(g1),gr_maxy(g1)
          write(*,*) '        minz1,maxz1 = ',gr_minz(g1),gr_maxz(g1)
          write(*,*) '        minx1,maxx2 = ',gr_minx(g2),gr_maxx(g2)
          write(*,*) '        miny1,maxy2 = ',gr_miny(g2),gr_maxy(g2)
          write(*,*) '        minz1,maxz2 = ',gr_minz(g2),gr_maxz(g2)
          write(*,*) '                        si, fi = ',si, fi
          write(*,*) '                        sj, fj = ',sj, fj
          write(*,*) '                        sk, fk = ',sk, fk
          write(*,*) '                        Length = ',length
          write(*,*) 'grid_find_intersection: Done.'
      end if

      return
      end       ! END: grid_find_intersection

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_return_owner:                                                        cc
cc                  Returns process which owns this grid.                     cc
cc                  NB: Now that the coarse level can have multiple grids,    cc
cc                      the coarse grid has a virtual parent owned by the     cc
cc                      master process.                                       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function grid_return_owner(grid)
      implicit    none
      include    'grid.inc'
      include    'mpif.h'
      include    'mpi_stuff.inc'
      logical     grid_return_existence
      external    grid_return_existence
      integer     grid

      !
      ! For example, the coarse grid has no parent
      ! so grid_return_owner(grid_return_parent(grid)) = -1
      !
      if (grid_return_existence(grid)) then
         grid_return_owner = gr_own(grid)
      else
         grid_return_owner = master
      end if

      return
      end       ! END: grid_return_owner

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_return_count:                                                        cc
cc                  Returns counter for grid.                                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function grid_return_count(gridnum)
      implicit    none
      integer     gridnum
      include    'grid.inc'

      grid_return_count = gr_count(gridnum)

      return
      end       ! END: grid_return_count

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_return_level:                                                        cc
cc                  Returns level of grid.                                    cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function grid_return_level(gridnum)
      implicit    none
      integer     gridnum
      include    'grid.inc'

      grid_return_level = gr_level(gridnum)

      return
      end       ! END: grid_return_level


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_return_parent:                                                       cc
cc                  Returns parent of grid.                                   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function grid_return_parent(gridnum)
      implicit    none
      integer     gridnum
      include    'grid.inc'

      grid_return_parent = gr_parent(gridnum)

      return
      end       ! END: grid_return_parent

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_pending:                                                             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function grid_pending(gridnum)
      implicit    none
      integer     gridnum
      include    'grid.inc'

      grid_pending = .false.
      if ( gr_alive(gridnum).eq. PENDING ) then
         grid_pending = .true.
      end if

      return
      end       ! END: grid_pending

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_dead:                                                                cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function grid_dead(gridnum)
      implicit    none
      integer     gridnum
      include    'grid.inc'

      grid_dead = .false.
      if ( gr_alive(gridnum).eq.DEAD ) then
         grid_dead = .true.
      end if

      return
      end       ! END: grid_dead

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_alive:                                                               cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function grid_alive(gridnum)
      implicit    none
      integer     gridnum
      include    'grid.inc'

      grid_alive = .false.
      if ( gr_alive(gridnum).eq.ALIVE ) then
         grid_alive = .true.
      end if

      return
      end       ! END: grid_alive

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_return_existence:                                                    cc
cc                  Returns true if grid exists, false if not.                cc
cc                  Does not actually test if grid has been created           cc
cc                  by searching through the hierarchy. Instead, it           cc
cc                  basically just differentiates a "null grid pointer"       cc
cc                  from a real (or possibly real) one.                       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function grid_return_existence(gridnum)
      implicit    none
      integer     gridnum
      include    'grid.inc'

      if ( (gridnum .gt. 0) .and. (gridnum .le. maxnumgrids) ) then
         grid_return_existence = .true.
      else
         grid_return_existence = .false.
      end if

      return
      end       ! END: grid_return_existence

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_return_child:                                                        cc
cc                  Returns pointer to first child of input grid.             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function grid_return_child(gridnum)
      implicit    none
      integer     gridnum
      include    'grid.inc'

      grid_return_child = gr_child(gridnum)

      return
      end       ! END: grid_return_child

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_return_grandchild:                                                   cc
cc                  Returns pointer to a grandchild of input grid.            cc
cc                  Used mostly to see if a grandchild exists.                cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function grid_return_grandchild(gridnum)
      implicit    none
      integer     gridnum
      include    'grid.inc'
      integer     gi_tmp
      logical     grid_return_existence
      logical     ltrace
      parameter ( ltrace = .false. )

      if (ltrace) write(*,*) 'grid_return_grandchild: Enter gi:',gridnum
      !
      ! Initialize to no grandchildren
      !
      grid_return_grandchild = -1

      gi_tmp = gr_child(gridnum)
      !
      ! Search for children of any of its children:
      !
 10   if (  grid_return_existence(gi_tmp)  .and.
     *     (gr_parent(gi_tmp).eq.gridnum )        ) then
         grid_return_grandchild = gr_child(gi_tmp)
         if (ltrace) write(*,*) ' g_r_g: ',grid_return_grandchild
         if (grid_return_existence(grid_return_grandchild)) goto 20
         gi_tmp = gr_sibling(gi_tmp)
         goto 10
      end if

 20   continue
      if (ltrace) write(*,*) 'grid_return_grandchild: Done. gc = ',
     *                        grid_return_grandchild

      return
      end       ! END: grid_return_grandchild

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_reset_time:                                                          cc
cc                  Sets grid time to zero.                                   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_reset_time(gridnum)
      implicit    none
      integer     gridnum
      include    'grid.inc'
      !include    'glob.inc'
      !integer     proc_return_myid

      logical     ltrace
      parameter ( ltrace = .false. )

      gr_t(gridnum) = 0.d0

      return
      end       ! END: grid_reset_time

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_advance_time:                                                        cc
cc                  Adds lambda*gr_h() to gr_time.                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_advance_time(gridnum)
      implicit    none
      integer     gridnum
      include    'grid.inc'
      include    'glob.inc'
      integer     proc_return_myid

      logical     ltrace
      parameter ( ltrace = .false. )

      gr_t(gridnum) = gr_t(gridnum) + lambda * gr_h(gridnum)

      if (ltrace) write(*,99) proc_return_myid(),
     *                        'grid_advance_time: old-->new: ',
     *                        gridnum,
     *                        gr_t(gridnum)-lambda*gr_h(gridnum),
     *                        gr_t(gridnum)
 99   format('[',I3,']',A,I5,2F12.4)

      return
      end       ! END: grid_advance_time

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_set_time:                                                            cc
cc                  Sets the time of a grid to that of its parent level.      cc
cc                  This is meant to fix any small (10e-14) disparaties       cc
cc                  when levels are time aligned (as determined by integer    cc
cc                  step counters)                                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_set_time(gridnum,time)
      implicit    none
      integer     gridnum
      real(kind=8) time  
      include    'grid.inc'
      include    'glob.inc'
      integer     proc_return_myid

      logical     ltrace
      parameter ( ltrace = .false. )

      gr_t(gridnum) = time

      return
      end       ! END: grid_set_time

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_return_minx:                                                         cc
cc                  Returns value of the x-anchor for this grid.              cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8) function grid_return_minx(gridnum)
      implicit    none
      integer     gridnum
      include    'grid.inc'

      grid_return_minx = gr_minx(gridnum)

      return
      end       ! END: grid_return_minx

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_return_miny:                                                         cc
cc                  Returns value of the y-anchor for this grid.              cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8) function grid_return_miny(gridnum)
      implicit    none
      integer     gridnum
      include    'grid.inc'

      grid_return_miny = gr_miny(gridnum)

      return
      end       ! END: grid_return_miny

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_return_minz:                                                         cc
cc                  Returns value of the z-anchor for this grid.              cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8) function grid_return_minz(gridnum)
      implicit    none
      integer     gridnum
      include    'grid.inc'

      grid_return_minz = gr_minz(gridnum)

      return
      end       ! END: grid_return_minz

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_return_time:                                                         cc
cc                  Returns time to which grid is advanced.                   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8) function grid_return_time(gridnum)
      implicit    none
      integer     gridnum
      include    'grid.inc'

      grid_return_time = gr_t(gridnum)

      return
      end       ! END: grid_return_time

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_dump_freenums:                                                       cc
cc                  Dump a listing of free grid numbers.                      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_dump_freenums()
      implicit    none
      include    'grid.inc'
      integer     proc_return_myid
      external    proc_return_myid
      integer     gi, myid

      myid = proc_return_myid()

      write(*,99) myid, '            numgrids = ',numgrids
      if (num_freegrids.eq.0) return

      write(*,*) '          //////\\\\\\'
      do gi = 1, num_freegrids
         write(*,99) myid,
     *           '   free grid number:  ',gi,free_grid_numbers(gi)
      end do
      write(*,*) '          \\\\\\//////'

  99  format('[',I3,'] ',A,2I6)

      return
      end       ! END: grid_dump_freenums

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_return_number_norec:                                                 cc
cc                  Returns a free grid number but without                    cc
cc                  reserving the number.                                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function grid_return_number_norec()
      implicit    none
      include    'grid.inc'
      integer     gridnum

      logical     ltrace
      parameter ( ltrace = .false. )


      !
      ! get "free" grid number for grid:
      !
      if (num_freegrids .ge. 1 ) then
         !  grab a previously used and then freed number:
         gridnum          = free_grid_numbers(num_freegrids)
         if (ltrace) then
            write(*,*) 'grid_return_number_norec: Using free'
            write(*,*) 'grid_return_number_norec: num_freegrids =',
     *                                            num_freegrids
         end if
      else
         ! take next grid number:
         if (numgrids .ge. maxnumgrids) then
            write(*,*) 'grid_return_number: No free grids available.'
            write(*,*) 'grid_return_number: Increase maxnumgrids,recomp'
            write(*,*) 'grid_return_number: maxnumgrids= ',maxnumgrids
            write(*,*) 'grid_return_number: Quitting...'
            call my_exit('No free grid#s available')
         end if
         gridnum          = numgrids + 1
         if (ltrace) write(*,*) 'grid_return_number_norec: getting next'
      end if

      if (ltrace) write(*,*)'grid_return_number_norec: Number: ',gridnum

      grid_return_number_norec = gridnum

      return
      end       ! END: grid_return_number_norec

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_return_number:                                                       cc
cc                  Returns a free grid number.                               cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function grid_return_number()
      implicit    none
      include    'grid.inc'
      integer     gridnum

      logical     ltrace
      parameter ( ltrace = .false. )

      !
      ! get "free" grid number for grid:
      !
      if (num_freegrids .ge. 1 ) then
         !  grab a previously used and then freed number:
         if (ltrace) then
            write(*,*) 'grid_return_number: Using free'
            write(*,*) 'grid_return_number: num_freegrids =',
     *                                            num_freegrids
         end if
         gridnum          = free_grid_numbers(num_freegrids)
         num_freegrids    = num_freegrids - 1
      else
         ! take next grid number:
         if (numgrids .ge. maxnumgrids) then
            write(*,*) 'grid_return_number: No free grids available.'
            write(*,*) 'grid_return_number: Increase maxnumgrids,recomp'
            write(*,*) 'grid_return_number: maxnumgrids= ',maxnumgrids
            write(*,*) 'grid_return_number: Quitting...'
            call my_exit('No free grid#s available')
         end if
         if (ltrace) write(*,*) 'grid_return_number: getting next'
         numgrids         = numgrids + 1
         gridnum          = numgrids
         if (numgrids .gt. max_numgrids) max_numgrids = numgrids
      end if

      grid_return_number = gridnum

      return
      end       ! END: grid_return_number

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_return_sibling:                                                      cc
cc                  Returns pointer to sibling of input grid.                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function grid_return_sibling(gridnum)
      implicit    none
      integer     gridnum
      include    'grid.inc'

      grid_return_sibling = gr_sibling(gridnum)

      return
      end       ! END: grid_return_sibling

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_set_owner:                                                           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_set_owner(grid, owner)
      implicit    none
      integer     grid,   owner
      include    'grid.inc'

      gr_own(grid) = owner

      return
      end          ! END: grid_set_owner

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_newA:                                                                cc
cc            NO anchors nor unique parents.                                  cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_newA(gi,own,lev,nx,ny,nz,minx,miny,minz)
      implicit      none
      integer       gi,     own,  lev, 
     *              nx,     ny,   nz
      real(kind=8)  minx,miny,minz
      include    'grid.inc'
      include    'glob.inc'
      include    'param.inc'
      include    'grid_methods.inc'
      integer     mem_alloc, gridnum, i, tmp_gi, myid
      real(kind=8)      parent_h, parent_minx, parent_miny, parent_minz
      integer     proc_return_myid
      external    proc_return_myid

      logical     ltrace
      parameter ( ltrace = .false. )


      myid = proc_return_myid()

      if (ltrace) then
         write(*,99)myid,'grid_newA: -------'
         write(*,99)myid,'grid_newA: Adding grid '
         write(*,99)myid,'grid_newA: gi      = ', gi
         write(*,99)myid,'grid_newA: own     = ', own
         write(*,99)myid,'grid_newA: lev     = ', lev
         write(*,99)myid,'grid_newA: nx      = ', nx
         write(*,99)myid,'grid_newA: ny      = ', ny
         write(*,99)myid,'grid_newA: nz      = ', nz
      end if
      if (ltrace) then
         write(*,99)myid,'grid_newA:ARun level_test on'
         call level_test_local()
         !call load_pointers(gi)
      end if

      !
      ! Store info about grid:
      !
      gr_level(gi)       = lev
      gr_own(gi)         = own
      gr_count(gi)       = 0
      gr_nx(gi)          = nx
      gr_ny(gi)          = ny
      gr_nz(gi)          = nz

      !
      ! Grids can be dead (and at some point removed from memory).
      ! During regridding process, grids can be "kill"ed but kept
      ! in memory to use their data before they're freed.
      ! Grid are created alive:
      !
      gr_alive(gi)     = ALIVE

      !
      ! This is for the RK3 routine...the grid needs to store
      ! which iteration it is, so we initialize it here:
      !
      gr_iter(gi)      = 0

      !
      ! Set the grid spacing:
      !
      if (lev .eq. 0) then
         gr_h(gi)      = h
      else
         gr_h(gi)      = grid_return_h(level_return_start(lev-1))
     *                    / (1.d0*refine_factor)
      end if

      !
      ! Set grid bounds:
      !
      gr_minx(gi)   = minx
      gr_miny(gi)   = miny
      gr_minz(gi)   = minz

      gr_maxx(gi)   = gr_minx(gi) + (nx-1)*gr_h(gi)
      gr_maxy(gi)   = gr_miny(gi) + (ny-1)*gr_h(gi)
      gr_maxz(gi)   = gr_minz(gi) + (nz-1)*gr_h(gi)

      !
      ! No children initially
      !
      gr_child(gi)  = -1

      !
      ! Set time of grid:
      !
      if (lev .eq. 0) then
         ! Must be the coarse grid, so we know time=0
         gr_t(gi)        = 0.d0
      else
         gr_t(gi)        = grid_return_time(level_return_start(lev-1))
      end if


         gr_work(gi) =    1.d0*( gr_nx(gi)    *gr_ny(gi)    *gr_nz(gi))
     *                  /(1.d0*( nx0          *ny0          *nz0      ))
         gr_work(gi) =    gr_work(gi)*refine_factor**lev

      if (ltrace) then
         write(*,97)myid,'grid_new: gr_work = ', gr_work(gi)
         write(*,97)myid,'grid_new: gr_t    = ', gr_t(gi)
         write(*,97)myid,'grid_new: gr_h    = ', gr_h(gi)
         write(*,97)myid,'grid_new: gr_minx = ', gr_minx(gi)
         write(*,97)myid,'grid_new: gr_miny = ', gr_miny(gi)
         write(*,97)myid,'grid_new: gr_minz = ', gr_minz(gi)
         write(*,97)myid,'grid_new: gr_maxx = ', gr_maxx(gi)
         write(*,97)myid,'grid_new: gr_maxy = ', gr_maxy(gi)
         write(*,97)myid,'grid_new: gr_maxz = ', gr_maxz(gi)
      end if

      !
      ! Update workload for owner process:
      !
      if(ltrace)write(*,99)myid,'grid_new: Calling proc_addwork'
      call proc_addwork(own,gr_work(gi), nx*ny*nz)

      !
      ! Need to update either "numgrids" or the list of free grid numbers
      !
         tmp_gi = grid_return_number()
         if (tmp_gi .ne. gi) then
            write(*,99)myid,'Problem Problem Problem'
            write(*,99)myid,'tmp_gi = ',tmp_gi
            write(*,99)myid,'gi     = ',gi
            write(*,99)myid,'own    = ',own
            write(*,99)myid,'myid   = ',myid
            write(*,99)myid,'When updating grid no,return value not exp'
            call level_tree_dump()
            call grid_dump_freenums()
            call my_exit('grid_newA:Updating grid number not the same')
         end if

      !
      ! Put grid into tree:
      !
      if (Levelp(lev).lt.0) then
         Levelp(lev)    = gi
         gr_sibling(gi) = -1
      else
         gr_sibling(gi) = Levelp(lev)
         Levelp(lev)    = gi
      end if

      if (ltrace) then
         !  Cannot run this test since the field pointers are not there yet:
         !write(*,99)myid,'grid_newA:SRun level_test on'
         !call level_test_local()
         !call load_pointers(gi)
         write(*,99)myid,'grid_newA: Done.'
      end if

 97   format('[',I4,'] ',A,3G15.5)
 99   format('[',I4,'] ',A,3I7)

      return
      end          ! END: grid_newA

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_createA:                                                             cc
cc                       Does not use anchor points.                          cc
cc                       Nor require a unique parent.                         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_createA(gi,own,lev,nx,ny,nz,minx,miny,minz)
      implicit    none
      integer      gi,     own,  lev,
     *             nx,     ny,   nz
      real(kind=8) minx,miny,minz
      include    'grid.inc'
      include    'glob.inc'
      !include    'param.inc'    ! for "use_mask"
      include    'mask.inc'
      integer     i, size
      integer     mem_alloc
      external    mem_alloc
      integer     proc_return_myid, myid
      external    proc_return_myid
      logical     ltrace
      parameter ( ltrace = .false. )

      myid = proc_return_myid()

      if (gi .le. 0) then
         write(*,99)myid,'grid_createA: Problem, gi must be > 0, gi:',gi
         call my_exit('grid_createA: Problem, gi must be > 0')
      end if

      if (ltrace) then
         write(*,99)myid,'grid_createA: Creating grid, gi: ',gi
         write(*,99)myid,'                gi    = ', gi
         write(*,99)myid,'                own   = ', own
         write(*,99)myid,'                lev   = ', lev
         write(*,99)myid,'                nx    = ', nx
         write(*,99)myid,'                ny    = ', ny
         write(*,99)myid,'                nz    = ', nz
         call mem_stat()
      end if

      !
      ! Establish grid
      !
      if (ltrace) then
         write(*,99)myid,'grid_createA: calling grid_newA'
         write(*,99)myid,'grid_createA:BRun level_test on'
         call level_test_local()
         !call load_pointers(gi)
      end if
      call grid_newA(gi, own, lev,nx,ny,nz,minx,miny,minz)

      !
      ! Allocate space and store pointers for fields:
      !
!     if (ltrace) then
!        write(*,99)myid,'grid_createA:QRun level_test on'
!        call level_test_local()
!        !call load_pointers(gi)
!     end if
      gr_x(gi)        = mem_alloc(nx + ny + nz)! NB:these 3 have to be
      if (gr_x(gi) .lt. 0) then
         write(*,99)myid,'grid_createA: Not enough memory for coords'
         return
      end if
      gr_y(gi)        = gr_x(gi) + nx      ! contiguous in memory
      gr_z(gi)        = gr_y(gi) + ny      ! used in gft_out_full call

!     if (ltrace) then
!        write(*,99)myid,'grid_createA:RRun level_test on'
!        call level_test_local()
!        !call load_pointers(gi)
!     end if
      size = nx*ny*nz
      do i = 1, num_gfuncs
         gfunc_pointer(i,gi)      = mem_alloc(size)
         if (gfunc_pointer(i,gi) .lt. 0) then
            write(*,99)myid,'grid_createA: Not enough memory field:',i
            return
         end if
         gfunc_tmpptr(i)          = gfunc_pointer(i,gi)
         if(ltrace) write(*,98)myid,'grid_createA: ptrs:',i,
     *          gfunc_pointer(i,gi),gfunc_name(i)
      end do

      if(ltrace)write(*,99)myid,'grid_createA: initializling chr'
      call load_pointers(gi)
      !
      ! Initialize CHR array:
      !
      !if (ltrace) then
         !write(*,99)myid,'grid_createA:SRun level_test on'
         !call level_test_local()
         !call load_pointers(gi)
      !end if
      call grid_init_chr(   q(gr_chr), 
     *                      gr_minx(gi), gr_miny(gi),
     *                      gr_minz(gi), gr_h(gi),
     *                      gr_nx(gi),   gr_ny(gi),
     *                      gr_nz(gi)                      )

      if(ltrace)write(*,99)myid,'grid_create: initializing coordinates'

      ! Initialize coordinates x(), y() and z():
      !if (ltrace) then
         !write(*,99)myid,'grid_createA:TRun level_test on'
         !call level_test_local()
         !call load_pointers(gi)
      !end if
      call linear_ramp( q(gr_x(gi)), gr_minx(gi),
     *                  gr_h(gi),    gr_nx(gi)    )
      call linear_ramp( q(gr_y(gi)), gr_miny(gi),
     *                  gr_h(gi),    gr_ny(gi)    )
      call linear_ramp( q(gr_z(gi)), gr_minz(gi),
     *                  gr_h(gi),    gr_nz(gi)    )


      if(ltrace)write(*,99)myid,'grid_create: initializing r array'

      ! Initialize "r" array to radial distance from origin:
      !if (ltrace) then
         !write(*,99)myid,'grid_createA:URun level_test on'
         !call level_test_local()
         !call load_pointers(gi)
      !end if
      call load_r(     q(gr_r), gr_h(gi),
     *        gr_minx(gi), gr_miny(gi), gr_minz(gi), nx,ny,nz )

      if(ltrace)write(*,99)myid,'grid_create: initializing flag array'

      ! Initialize "flag" array to zero:
      !if (ltrace) then
         !write(*,99)myid,'grid_createA:VRun level_test on'
         !call level_test_local()
         !call load_pointers(gi)
      !end if
      !call load_scal3d(q(gr_flag), FLAG_NOREFINE, nx,ny,nz )
      call load_scal3d(q(gr_flag), 0.d0, nx,ny,nz )
      call load_scal3d(q(gr_tmp),  0.d0, nx,ny,nz )

      ! Initialize to all unmasked points:
      !if (ltrace) then
         !write(*,99)myid,'grid_createA:YRun level_test on'
         !call level_test_local()
         !call load_pointers(gi)
      !end if
      call load_scal3D(q(gr_mask),1.d0*BHMASK_UNMASKED,nx,ny,nz)
      call load_scal3D(q(gr_wdiss),1.d0,nx,ny,nz)

      !if (ltrace) then
         !write(*,99)myid,'grid_createA:ZRun level_test on'
         !call level_test_local()
      !end if
      if(ltrace)write(*,99)myid,'grid_createA: Done creating grid'

 98   format('[',I4,'] ',A,2I10,A)
 99   format('[',I4,'] ',A,3I7)
      return
      end     ! END: subroutine grid_createa


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_inject_master:                                                       cc
cc                 For MPI version for remotely injecting fields.             cc
cc                 Receives correction for processors which own children      cc
cc                 grids and applies those corrections to the respective      cc
cc                 fields. Complements grid_inject_slave.                     cc
cc             NB: Updated for overlapping grid using CHR mask for averaging. cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_inject_master(child,parent)
      implicit    none
      integer     child, parent
      include    'grid.inc'
      include    'glob.inc'
      include    'mpif.h'
      include    'mpi_stuff.inc'
      integer     length, owner
      integer     nxp, nyp, nzp, nxc, nyc, nzc, i, nxs,nys,nzs
      integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
      integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
      integer     tmp_ghostwidth
      logical     double_equal
      external    double_equal
      character(3) tmps,tmps2
      ! Field number to be output if tracing is on:
      integer     outfnum
      parameter ( outfnum  = 1 )
      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )

      !ltrace  =.false.
      !ltrace2 =.false.
      !if (child.eq.19) then
      !ltrace  =.true.
      !ltrace2 =.true.
      !end if

      !call grid_check_mask(parent,20)

      owner  = gr_own(child)

      !
      ! This is the grid this process owns
      !
      call load_pointers(parent)

      !
      ! Parent's dimensions:
      !
      nxp = gr_nx(parent)
      nyp = gr_ny(parent)
      nzp = gr_nz(parent)

      !
      ! Child's dimensions:
      !
      nxc = gr_nx(child)
      nyc = gr_ny(child)
      nzc = gr_nz(child)

      !
      ! After an elliptic solve initially, we inject but want
      ! to inject with ghostwidth==1 no matter what the real
      ! value is:
      !
      if (double_equal(gr_t(child),0.d0)) then
         tmp_ghostwidth = 1
      else
         tmp_ghostwidth = ghostwidth
      end if

      call grid_find_intersection(child, parent, length,
     *            min_ic, max_ic, min_jc, max_jc, min_kc, max_kc)
      call grid_find_intersection(parent, child, length,
     *            min_ip, max_ip, min_jp, max_jp, min_kp, max_kp)

      !
      ! Overlap region on the parent
      !    but with resolution of the parent grid:
      !
      nxs = max_ip - min_ip + 1
      nys = max_jp - min_jp + 1
      nzs = max_kp - min_kp + 1

      !
      ! (linear) Size of grid to be transmitted:
      !
      length = nxs*nys*nzs

      if (ltrace) then
        call int2str(child,tmps)
        call int2str(parent,tmps2)
         write(*,99) myid,'grid_inject_master: child   = ', child
         write(*,99) myid,'grid_inject_master: parent  = ', parent
         write(*,99) myid,'grid_inject_master: length  = ', length
         write(*,99) myid,'grid_inject_master: owner   = ', owner
         write(*,99) myid,'grid_inject_master: nx/y/zp = ',nxp,nyp,nzp
         write(*,99) myid,'grid_inject_master: nx/y/zp = ',nxs,nys,nzs
         write(*,99) myid,'grid_inject_master: nx/y/zc = ',nxc,nyc,nzc
         write(*,99) myid,'grid_inject_master: mini/j/kp=',
     *                                            min_ip,min_jp,min_kp
         write(*,99) myid,'grid_inject_master: mini/j/kc=',
     *                                            min_ic,min_jc,min_kc
         write(*,99) myid,'grid_inject_master: maxi/j/kc=',
     *                                            max_ic,max_jc,max_kc
      end if

      call load_scal3d(q(gr_tmp), 0.d0, nxp,nyp,nzp )

      do i = 1, num_gfuncs
         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL .or.
     *        gfunc_name(i) .eq. 'mask'         .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL.or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRAL_DER ) then
            !
            ! Receive correction for grid owner:
            !
            if (ltrace2) then
               write(*,99) myid,' master: receiving correction'
               write(*,99) myid,' master:     field: ',i
               write(*,99) myid,' master: from owner: ', owner
            end if
            call MPI_Recv(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_INJ_DATA, MPI_COMM_WORLD, status,ierr)
            if (i.eq.outfnum .and. ltrace)
     *      call field_out3d(q(gr_tmp),gr_t(parent),
     *           'InjectMTmp'//tmps//tmps2,
     *       q(gr_x(parent)+(min_ip-1)),q(gr_x(parent)+(max_ip-1)),
     *       q(gr_y(parent)+(min_jp-1)),q(gr_y(parent)+(max_jp-1)),
     *       q(gr_z(parent)+(min_kp-1)),q(gr_z(parent)+(max_kp-1)),
     *             nxs,nys,nzs,myid)
            !
            ! Replace values from tmp storage:
            !
            call replace_part_field_avgGWc( 
     *                           q(gfunc_pointer(i,parent)),
     *                           q(gr_chr),
     *                           q(gr_tmp),
     *                           min_ip,
     *                           min_jp,
     *                           min_kp,
     *                           nxp,nyp,nzp,
     *                           nxs,nys,nzs,
     *                           nxc,nyc,nzc,
     *                           min_ic,max_ic,
     *                           min_jc,max_jc,
     *                           min_kc,max_kc,
     *                           refine_factor, tmp_ghostwidth,
     *                           bound_width )
            if (ltrace2)
     *      call field_dump_info(q(gfunc_pointer(i,parent)),nxs,nys,nzs)
            if (i.eq.outfnum .and. ltrace)
     *      call field_out3d(q(gfunc_pointer(i,parent)),gr_t(parent),
     *                   'InjectMPar'//tmps//tmps2,
     *               gr_minx(parent),gr_maxx(parent),
     *               gr_miny(parent),gr_maxy(parent),
     *               gr_minz(parent),gr_maxz(parent),
     *               nxp,nyp,nzp,myid)
         else if ( gfunc_type(i) .eq. GFUNC_INTEGRALFV .or.
     *             gfunc_type(i) .eq. GFUNC_INTEGRAL_DERFV ) then
            call MPI_Recv(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_INJ_DATA, MPI_COMM_WORLD, status,ierr)
            !call load_scal1d(q(gr_tmp),-1.d0*child,nxs*nys*nzs)
            call replace_part_field_avgFV(
     *                           q(gfunc_pointer(i,parent)),
     *                           q(gr_chr),
     *                           q(gr_tmp),
     *                           min_ip,
     *                           min_jp,
     *                           min_kp,
     *                           nxp,nyp,nzp,
     *                           nxs,nys,nzs,
     *                           nxc,nyc,nzc,
     *                           min_ic,max_ic,
     *                           min_jc,max_jc,
     *                           min_kc,max_kc,
     *                           refine_factor, tmp_ghostwidth,
     *                           bound_width )
         end if
      end do

      !
      ! this isn't field specific, but we can use the last field
      ! updated above to get info about which points are being
      ! restricted:
      !
      call inc_mask_overlapB(    q(gr_chr),
     *                           q(gr_tmp),
     *                           min_ip,
     *                           min_jp,
     *                           min_kp,
     *                           nxp,nyp,nzp,
     *                           nxs,nys,nzs,
     *                           nxc,nyc,nzc,
     *                           min_ic,max_ic,
     *                           min_jc,max_jc,
     *                           min_kc,max_kc,
     *                           refine_factor, tmp_ghostwidth )

      ! Enforce symmetry conditions after injection:
      if (assume_symmetry.ne.0) call grid_symm(parent)

      if (ltrace)
     *      call field_out3d(q(gfunc_pointer(outfnum,parent)),
     *                 gr_t(parent),
     *                   'InjectMPOst'//tmps//tmps2,
     *               gr_minx(parent),gr_maxx(parent),
     *               gr_miny(parent),gr_maxy(parent),
     *               gr_minz(parent),gr_maxz(parent),
     *               nxp,nyp,nzp,myid)

      !call grid_check_mask(parent,21)

      call load_scal3d(q(gr_tmp), 0.d0, nxp,nyp,nzp )

      if (ltrace) then
         call grid_dump_info(parent)
         !call grid_dump_info(child)
         write(*,99)myid,'grid_inject_master: Done.',parent, child
      end if
 99   format('[',I3,'] ',A,3I7)

      return
      end    ! END: grid_inject_master

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_inject_slave:                                                        cc
cc                 For MPI version for remotely injecting fields.             cc
cc                 Restricts fields into tmp storage, and then passes that    cc
cc                 to the parent of this grid so that "grid_inject_master"    cc
cc                 does the rest of the work in terms of updating the fields. cc
cc             NB: Apparently no change to allow for overlapping grids.       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_inject_slave(child,parent)
      implicit    none
      integer     child, parent
      include    'grid.inc'
      include    'glob.inc'
      include    'mpif.h'
      include    'mpi_stuff.inc'
      integer     length, owner
      integer     nxp, nyp, nzp, nxc, nyc, nzc, i, nxs,nys,nzs
      integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
      integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
      integer     tmp_ghostwidth
      character(3) tmps,tmps2
      logical     double_equal
      external    double_equal
      ! Field number to be output if tracing is on:
      integer     outfnum
      parameter ( outfnum  = 1 )
      logical     ltrace
      !parameter ( ltrace  = .false. )
      logical     ltrace2
      !parameter ( ltrace2 = .false. )

      ltrace  =.false.
      ltrace2 =.false.
!     if (child.eq.179) then
!     ltrace  =.true.
!     ltrace2 =.true.
!     end if

      owner  = gr_own(parent)

      call load_pointers(child)

      !
      ! After an elliptic solve initially, we inject but want
      ! to inject with ghostwidth==1 no matter what the real
      ! value is:
      !
      if (double_equal(gr_t(child),0.d0)) then
         tmp_ghostwidth = 1
      else
         tmp_ghostwidth = ghostwidth
      end if

      !
      ! Parent's dimensions:
      !
      nxp = gr_nx(parent)
      nyp = gr_ny(parent)
      nzp = gr_nz(parent)

      !
      ! Child's dimensions:
      !
      nxc = gr_nx(child)
      nyc = gr_ny(child)
      nzc = gr_nz(child)

      call grid_find_intersection(child, parent, length,
     *            min_ic, max_ic, min_jc, max_jc, min_kc, max_kc)
      call grid_find_intersection(parent, child, length,
     *            min_ip, max_ip, min_jp, max_jp, min_kp, max_kp)

      !
      ! Overlap region on the parent
      !    but with resolution of the parent grid:
      !
      nxs = max_ip - min_ip + 1
      nys = max_jp - min_jp + 1
      nzs = max_kp - min_kp + 1

      !
      ! (linear) Size of grid to be transmitted:
      !
      length = nxs*nys*nzs

      if (ltrace) then
        call int2str(child,tmps)
        call int2str(parent,tmps2)
         write(*,91) myid,'grid_inject_slave: c/p/t= ', child,parent,
     .                                                    gr_t(child)
         write(*,99) myid,'grid_inject_slave: child   = ', child
         write(*,99) myid,'grid_inject_slave: parent  = ', parent
         write(*,99) myid,'grid_inject_slave: length  = ', length
         write(*,99) myid,'grid_inject_slave: owner   = ', owner
         write(*,99) myid,'grid_inject_slave: nx/y/zp = ',nxp,nyp,nzp
         write(*,99) myid,'grid_inject_slave: nx/y/zs = ',nxs,nys,nzs
         write(*,99) myid,'grid_inject_slave: nx/y/zc = ',nxc,nyc,nzc
         write(*,99) myid,'grid_inject_slave: mini/j/kp=',
     *                                            min_ip,min_jp,min_kp
         write(*,99) myid,'grid_inject_slave: mini/j/kc=',
     *                                            min_ic,min_jc,min_kc
         write(*,99) myid,'grid_inject_slave: maxi/j/kc=',
     *                                            max_ic,max_jc,max_kc
      end if

      call load_scal3d(q(gr_tmp), 0.d0, nxc,nyc,nzc )

      do i = 1, num_gfuncs
         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL .or.
     *        gfunc_name(i) .eq. 'mask'         .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL.or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRAL_DER ) then
            !
            ! Restrict to tmp storage:
            !
            if (ltrace2)call field_dump_info(q(gfunc_pointer(i,child)),
     *                                         nxc,nyc,nzc)
            if (i.eq.outfnum .and. ltrace)
     *      call field_out3d(q(gfunc_pointer(i,child)),gr_t(child),
     *                   'InjectSChi'//tmps//tmps2,
     *               gr_minx(child),gr_maxx(child),
     *               gr_miny(child),gr_maxy(child),
     *               gr_minz(child),gr_maxz(child),
     *               nxc,nyc,nzc,myid)
            if (i.eq.outfnum .and. ltrace)
     *      call field_out3d(q(gr_chr), gr_t(child),
     *                   'InjectSChr'//tmps//tmps2,
     *               gr_minx(child),gr_maxx(child),
     *               gr_miny(child),gr_maxy(child),
     *               gr_minz(child),gr_maxz(child),
     *               nxc,nyc,nzc,myid)
            call restrict_fieldGWb(q(gr_tmp),
     *                           q(gfunc_pointer(i,child)),
     *                           q(gr_chr),
     *                           nxs,nys,nzs,
     *                           nxc,nyc,nzc,
     *                           min_ic,max_ic,
     *                           min_jc,max_jc,
     *                           min_kc,max_kc,
     *                           refine_factor,tmp_ghostwidth)
            if (ltrace2)call field_dump_info(q(gr_tmp),nxs,nys,nzs)
            if (i.eq.outfnum .and. ltrace)
     *      call field_out3d(q(gr_tmp),gr_t(child),
     *           'InjectSTmp'//tmps//tmps2,
     *       q(gr_x(child)+(min_ic-1)),q(gr_x(child)+(max_ic-1)),
     *       q(gr_y(child)+(min_jc-1)),q(gr_y(child)+(max_jc-1)),
     *       q(gr_z(child)+(min_kc-1)),q(gr_z(child)+(max_kc-1)),
     *             nxs,nys,nzs,myid)
            !
            ! Send to owner of parent:
            !
            if (ltrace2) then
               write(*,99) myid,' slave:  sending correction'
               write(*,99) myid,' slave:      field: ',i
               write(*,99) myid,' slave:   to owner: ', owner
            end if
            call MPI_Send(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_INJ_DATA, MPI_COMM_WORLD, ierr)
         else if ( gfunc_type(i) .eq. GFUNC_INTEGRALFV .or.
     *             gfunc_type(i) .eq. GFUNC_INTEGRAL_DERFV ) then
            call restrict_fieldFV(q(gr_tmp),
     *                           q(gfunc_pointer(i,child)),
     *                           q(gr_chr),
     *                           nxs,nys,nzs,
     *                           nxc,nyc,nzc,
     *                           min_ic,max_ic,
     *                           min_jc,max_jc,
     *                           min_kc,max_kc,
     *                           refine_factor,tmp_ghostwidth)
            call MPI_Send(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_INJ_DATA, MPI_COMM_WORLD, ierr)
         end if
      end do

      call load_scal3d(q(gr_tmp), 0.d0, nxc,nyc,nzc )

      if (ltrace) then
         !call grid_dump_info(parent)
         call grid_dump_info(child)
         !call grid_test(child)
         call level_test_local()
         write(*,99)myid,'grid_inject_slave: Done.',parent, child
      end if
 91   format('[',I3,'] ',A,2I7,G14.7)
 99   format('[',I3,'] ',A,3I7)

      return
      end    ! END: grid_inject_slave

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_reset_mask:                                                          cc
cc                 Reset the mask on a given grid.                            cc
cc                 The mask is used when child grids are injected into        cc
cc                 parents so that the injection can be averaged among        cc
cc                 the various children. Hence, the mask needs to be reset    cc
cc                 after injection.                                           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_reset_mask(grid)
      implicit    none
      integer     grid
      include     'grid.inc'
      logical     ltrace
      parameter ( ltrace = .false. )

      call load_pointers(grid)

      call remove_negs(     q(gr_chr), 
     *                      gr_nx(grid),   gr_ny(grid),
     *                      gr_nz(grid)                      )


      return
      end    ! END: grid_reset_mask

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_syncbnd_slaveLG:                                                     cc
cc                 Sync the boundaries of a pair of grids using MPI.          cc
cc                 This routine is always called such that this process       cc
cc                 owns grid "gi1" but does *not* own grid "gi2".             cc
cc                 Also, this routine is always called on some other          cc
cc                 with arguments reversed (obviously) at the same time       cc
cc                 so there is no grid_syncbnd_master.                        cc
cc             NB: This version derives from grid_syncbndslave() but          cc
cc                 collects data for each field into a single buffer          cc
cc                 and sends that LARGE message. Matt A. tells me that        cc
cc                 that this will be faster than a bunch of smaller msgs.     cc
cc                 Not clear whether I need to allocate a buffer or if        cc
cc                 gr_tmp will be big enough.                                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_syncbnd_slaveLG(gi1, gi2)
      implicit    none
      integer     gi1, gi2
      include     'grid.inc'
      include     'grid_methods.inc'
      include     'glob.inc'
      include     'chr.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'
      integer     parent, nxs,nys,nzs, length
      integer     lx, rx, ly, ry, lz, rz
      integer     nx1, ny1, nz1, nx2, ny2, nz2, i,j,k,
     *            min_i1, max_i1,
     *            min_j1, max_j1,
     *            min_k1, max_k1,
     *            min_i2, max_i2,
     *            min_j2, max_j2,
     *            min_k2, max_k2,
     *            properlength, maxlength,
     *            length1, length2, owner, index
      integer     request_recv, request_send
      logical     fvflag
      real(kind=8) oldnorm1,newnorm1, nrmdiff1, time
      logical     double_equal
      external    double_equal
      integer     calc_bnd_lenb
      external    calc_bnd_lenb
      real(kind=8) myl2norm, myl2norm3d, mynorm
      external     myl2norm, myl2norm3d
      character(3) tmps,tmps2
      integer     outbuf, inbuf, sizebuf, fieldcount, offset
      integer     gr_chr_gi1
      ! Field number to be output if tracing is on:
      integer     outfnum
      parameter ( outfnum  = 2 )
      logical     buffertoosmall, allocatebuffers
      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )
      ! Trace basically just the MPI stuff:
      logical     ltraces
      parameter ( ltraces = .false. )
      integer      mem_alloc
      external     mem_alloc
      logical      proc_havemem
      external     proc_havemem

 98   format('[',I4,'] ',A,3G15.8)
 99   format('[',I4,'] ',A,4I15)

      if(ltraces)write(*,99)myid,'grid_syncbnd_slaveLG: gi/gj:', gi1,gi2
      if (ltrace) then
         write(*,99) myid,'grid_syncbnd_slaveLG: gi, gj: ', gi1,gi2
         call grid_dump_info(gi1)
         call grid_dump_info(gi2)
      end if

      !
      ! Get pointers to temporary array:
      !
      call load_pointers(gi1)
      gr_chr_gi1 = gr_chr
      owner = grid_return_owner(gi2)

      time = gr_t(gi1)

      nx1 = gr_nx(gi1)
      ny1 = gr_ny(gi1)
      nz1 = gr_nz(gi1)

      nx2 = gr_nx(gi2)
      ny2 = gr_ny(gi2)
      nz2 = gr_nz(gi2)

      call grid_find_intersection(gi1, gi2,length1,
     *            min_i1, max_i1, min_j1, max_j1, min_k1, max_k1)

      !
      ! If no intersection:
      !
      if (length1 .eq. 0) return
      !
      ! If grids do not overlap at their ends
      !         more than the bound_width:
      !
      if (    ((max_i1-min_i1).lt.  bound_width-1 .and.
     *                (max_i1.eq.nx1.or.min_i1.eq.1))
     *    .or.((max_j1-min_j1).lt.  bound_width-1 .and.
     *                (max_j1.eq.ny1.or.min_j1.eq.1))
     *    .or.((max_k1-min_k1).lt.  bound_width-1 .and.
     *                (max_k1.eq.nz1.or.min_k1.eq.1)) ) then
         if (ltrace) then
               write(*,99)myid,'grid_syncbnd_slave:not syncing',gi1,gi2
               call grid_dump_info(gi1)
               call grid_dump_info(gi2)
               write(*,99)myid,'grid_syncbnd_slave:outputting chrgi1'
               call int2str(gi1,tmps)
               call int2str(gi2,tmps2)
             call field_out3d(q(gr_chr),time,     'chrgi1'//tmps//tmps2,
     *             gr_minx(gi1),gr_maxx(gi1),
     *             gr_miny(gi1),gr_maxy(gi1),
     *             gr_minz(gi1),gr_maxz(gi1),
     *             nx1,ny1,nz1,myid)
            write(*,99)myid,'grid_syncbnd_slave: min_i,j,k1= ',min_i1,
     *                      min_j1,min_k1
            write(*,99)myid,'grid_syncbnd_slave: max_i,j,k1= ',max_i1,
     *                      max_j1,max_k1
         end if
         return
      end if

      call grid_find_intersection(gi2, gi1,length2,
     *            min_i2, max_i2, min_j2, max_j2, min_k2, max_k2)


      !
      ! Not needed anymore:
      lx = 1
      rx = 1
      ly = 1
      ry = 1
      lz = 1
      rz = 1

      !
      ! Above lengths don't work for just the boundary:
      !
      length1 = calc_bnd_lenb(min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *                        min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *                        lx,rx, ly,ry, lz,rz,
     *                        nx1,ny1,nz1, nx2,ny2,nz2,bound_width)
      length2 = calc_bnd_lenb(min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *                        min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *                        lx,rx, ly,ry, lz,rz,
     *                        nx2,ny2,nz2, nx1,ny1,nz1,bound_width)


      maxlength = length1
      if (length2 .gt. maxlength) maxlength = length2
      !
      ! If no "effective" intersection:
      !
      if (maxlength .eq. 0) return
      !
      ! If overlap is just one point wide, do not sync
      ! since those points can't be interior
      !
      if (  (min_i1.eq.max_i1 .and .max_i1.eq.nx1) .or.
     *      (min_j1.eq.max_j1 .and .max_j1.eq.ny1) .or.
     *      (min_k1.eq.max_k1 .and .max_k1.eq.nz1) .or.
     *      (min_i2.eq.max_i2 .and .max_i2.eq.nx2) .or.
     *      (min_j2.eq.max_j2 .and .max_j2.eq.ny2) .or.
     *      (min_k2.eq.max_k2 .and .max_k2.eq.nz2) ) then
        if(ltrace)write(*,99)myid,'grid_syncbnd_slLG:Notsyncing',gi1,gi2
        if(ltrace)call grid_dump_info(gi1)
        if(ltrace)call grid_dump_info(gi2)
        return
      end if

      !
      ! check size of buffer:
      !
      buffertoosmall  = .false.
      allocatebuffers = .false.
      sizebuf         = maxlength*num_synfields
      inbuf           = gr_tmp
      outbuf          = gr_flag
      if (sizebuf .gt. nx1*ny1*nz1) then
         if(ltrace2) then
         write(*,99)myid,'grid_syncbndLG: maxlength    =',maxlength
         write(*,99)myid,'grid_syncbndLG: num_synfields=',num_synfields
         write(*,99)myid,'grid_syncbndLG: size reqd    =',sizebuf
         write(*,99)myid,'grid_syncbndLG: size avail.1 =',nx1*ny1*nz1
         end if
         buffertoosmall = .true.
      end if
      !
      if (sizebuf .gt. nx2*ny2*nz2) then
         if(ltrace2) then
         write(*,99)myid,'grid_syncbndLG: maxlength    =',maxlength
         write(*,99)myid,'grid_syncbndLG: num_synfields=',num_synfields
         write(*,99)myid,'grid_syncbndLG: size reqd    =',sizebuf
         write(*,99)myid,'grid_syncbndLG: size avail.2 =',nx2*ny2*nz2
         end if
         buffertoosmall = .true.
      end if
      !
!     write(*,*)myid,proc_havemem(myid,maxlength*num_synfields),
!    *               proc_havemem(owner,maxlength*num_synfields)
      if ( buffertoosmall
     .              .and. 
     .          proc_havemem(myid,sizebuf)
     .              .and. 
     .          proc_havemem(owner,sizebuf) ) then
          if(ltrace) then
             call mem_stat()
             write(*,99)myid,'grid_syncbndLG: Allocating buffers'
          end if
          inbuf           = mem_alloc(sizebuf)
          outbuf          = mem_alloc(sizebuf)
          allocatebuffers = .true.
          buffertoosmall  = .false.
      else if (buffertoosmall) then
          if(ltrace) then
            write(*,99)myid,'grid_syncbndLG: Not enough memory to'
            write(*,99)myid,'grid_syncbndLG: allocate buffers'
          end if
          call mem_stat()
      end if
      if (buffertoosmall) then
         call mem_stat()
         write(*,99)myid,'grid_syncbndLG: Buffer is too small!!'
         write(*,99)myid,'grid_syncbndLG: calling grid_syncbnd_slave'
         call grid_syncbnd_slave(gi1,gi2)
      end if

      call load_scal1D(q(inbuf),  0.d0, sizebuf )
      call load_scal1D(q(outbuf), 0.d0, sizebuf )

      !
      ! Non-blocking recv first:
      !
      if(ltraces)write(*,99)myid,'  LG: doing IRecv from ',owner
      call MPI_IRecv(q(inbuf), sizebuf, MPI_DOUBLE_PRECISION,
     *           owner, TAG_SYNC_DATA, MPI_COMM_WORLD,request_recv,ierr)
      if(ltraces)write(*,99)myid,'  LG: IRecv   ierr=',ierr,request_recv

      !
      ! Construct out-going buffer with necessary data:
      !
      fieldcount = 0
      do i = 1, num_gfuncs
         !
         if ((gfunc_type(i) .eq.GFUNC_INTEGRAL_ADV)   .or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRAL  .and.(time.eq.0)).or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRAL_DER)
     *                 .or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRAL_DEREL)
     *                 .or.
     *       (gfunc_class(i) .eq. GFUNC_W)
     *                 .or.
     *       (gfunc_name(i) .eq. 'mask')
     *                 .or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRAL_ADVFV).or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRALFV.and.time.eq.0).or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRAL_DERFV) ) then
            !
            if ( gfunc_type(i) .eq. GFUNC_INTEGRALFV .or.
     *           gfunc_type(i) .eq. GFUNC_INTEGRAL_DERFV ) then
               ! Finite Volume fields get sync'ed differently:
                fvflag = .true.
            else
                fvflag = .false.
            end if
            !
            if (fieldcount.gt.num_synfields-1) then
               write(*,99)myid,'syncbnd_slaveLG: Buffer not big enough!'
               write(*,99)myid,'syncbnd_slaveLG:         i=',i
               write(*,99)myid,'syncbnd_slaveLG:num_gfuncs=',num_gfuncs
               write(*,99)myid,'syncbnd_slaveLG:fieldcount=',fieldcount
            write(*,99)myid,'syncbnd_slaveLGnum_synfields',num_synfields
               write(*,99)myid,'syncbnd_slaveLG: Quitting...'
               call my_exit('fieldcount bigger than num_synfields')
            end if
            !
            ! Copy data into tmp now that it's safe:
            !
            offset = maxlength*fieldcount
            if(ltraces)then
               write(*,99)myid,'LG:  Copying data from i =',i
               write(*,99)myid,'LG:           fieldcount = ',fieldcount
               write(*,99)myid,'LG:           offset     = ',offset
               write(*,99)myid,'LG:           sizebuf    = ',sizebuf
               write(*,99)myid,'LG:           maxlength  = ',maxlength
            end if
            if (offset .gt. sizebuf-maxlength) then
               write(*,99)myid,'syncbnd_slaveLG: Offset very wrong'
               write(*,99)myid,'LG:           fieldcount = ',fieldcount
               write(*,99)myid,'LG:           offset     = ',offset
               write(*,99)myid,'LG:           sizebuf    = ',sizebuf
               write(*,99)myid,'LG:           maxlength  = ',maxlength
            end if
            call copy_int_bndd(   q(outbuf+offset),
     *                           q(gfunc_pointer(i,gi1)),
     *                           q(gr_chr_gi1),
     *            min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *            min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *                           lx,rx, ly,ry, lz,rz,
     *                           length1, properlength,
     *            nx1,ny1,nz1, nx2,ny2,nz2, bound_width,fvflag)
            !
            fieldcount = fieldcount + 1
            !
         else
            if(ltrace2)write(*,99)myid,'syncbnd_slave:ignorefield i= ',i
         end if
         !
      end do

      !
      ! Send outgoing buffer now:
      !
      if(ltraces)write(*,99)myid,'  Done copying data to outbuf'
      if(ltraces)write(*,98)myid,'  Doing ISend', 
     .                         myl2norm(q(outbuf),sizebuf)
      call MPI_ISend(q(outbuf),sizebuf,MPI_DOUBLE_PRECISION,
     *           owner, TAG_SYNC_DATA,MPI_COMM_WORLD, request_send,ierr)
      if (ltrace2.or.ltraces) then
         write(*,99)myid,'  ISend ierr: ',ierr
         write(*,99)myid,'  ISend requestsd:',request_send
      end if

      !
      ! Wait for receive to complete before pasting data into fields:
      !
      if(ltraces)write(*,99)myid,' Waiting.. recv',request_recv
      call MPI_Wait(request_recv, status, ierr)
      if(ltraces)write(*,99)myid,'  Done waiting recv, ierr: ',ierr

      !
      ! Paste data from receive buffer into appropriate fields:
      !
      fieldcount = 0
      do i = 1, num_gfuncs
         !
         if ((gfunc_type(i) .eq.GFUNC_INTEGRAL_ADV)   .or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRAL  .and.(time.eq.0)).or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRAL_DER)
     *                 .or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRAL_DEREL)
     *                 .or.
     *       (gfunc_class(i) .eq. GFUNC_W)
     *                 .or.
     *       (gfunc_name(i) .eq. 'mask')
     *                 .or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRAL_ADVFV).or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRALFV.and.time.eq.0).or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRAL_DERFV) ) then
            !
            if(ltraces)write(*,99)myid,'  pasting field i=',i,fieldcount
            !
            if ( gfunc_type(i) .eq. GFUNC_INTEGRALFV .or.
     *           gfunc_type(i) .eq. GFUNC_INTEGRAL_DERFV ) then
               ! Finite Volume fields get sync'ed differently:
                fvflag = .true.
            else
                fvflag = .false.
            end if
            !
            if (ltraces) write(*,99) myid,'  Pasting data from ',owner
            offset = maxlength*fieldcount
            if (offset .gt. sizebuf-maxlength) then
               write(*,99)myid,'syncbnd_slaveLG: Offset very wrong'
            end if
            call paste_int_bndd(  q(inbuf+offset),
     *                           q(gfunc_pointer(i,gi1)),
     *                           q(gr_chr_gi1),
     *            min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *            min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *                           lx,rx, ly,ry, lz,rz,
     *                           nx2,ny2,nz2,
     *                           properlength,
     *                           nx1,ny1,nz1,bound_width,fvflag )
            !
            fieldcount = fieldcount + 1
            !
         else
            if(ltrace2)write(*,99)myid,'syncbnd_slaveLG:ignorefiel i=',i
         end if
         !
      end do

      !
      ! Wait for send to complete before finishing up and returning:
      !
      if(ltraces)write(*,99)myid,'Waiting..send',request_send
      call MPI_Wait(request_send, status, ierr)
      if(ltraces)write(*,99)myid,'  Wait on send done.',ierr

      if (ltrace2) then
       newnorm1 = myl2norm3d(q(gfunc_pointer(outfnum,gi1)),nx1,ny1,nz1)
       nrmdiff1 = abs(newnorm1-oldnorm1)/oldnorm1
       write(*,98)myid,'grid_syncbnd_slave: newnorm1  = ',newnorm1
       write(*,98)myid,'grid_syncbnd_slave: nrmdiff1  = ',nrmdiff1
       if (nrmdiff1 .gt. 0.1) then
          write(*,98)myid,'grid_syncbnd_slave: PROBLEM!!!!! '
          write(*,98)myid,'grid_syncbnd_slave:norms1:',newnorm1,oldnorm1
       end if
       call field_out3d(q(gfunc_pointer(outfnum,gi1)),time,
     *             'postslave'//tmps//tmps2,
     *             gr_minx(gi1),gr_maxx(gi1),
     *             gr_miny(gi1),gr_maxy(gi1),
     *             gr_minz(gi1),gr_maxz(gi1),
     *             nx1,ny1,nz1,myid)
      end if

      !
      ! If we needed to allocate memory, free it now:
      !
      if (allocatebuffers) then
         if(ltrace)then
            write(*,99)myid,'grid_syncbnd_slaveLG: Deallocate buffers'
         end if
         call mem_dealloc(inbuf,  sizebuf )
         call mem_dealloc(outbuf, sizebuf )
      end if

      if(ltrace)then
         write(*,99)myid,'grid_syncbnd_slaveLG:level_test_local',gi1,gi2
         call level_test_local()
         write(*,99)myid,'grid_syncbnd_slaveLG: Done:',gi1,gi2
      end if

      return
      end    ! END: grid_syncbnd_slaveLG

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_syncbnd_slave:                                                       cc
cc                 Sync the boundaries of a pair of grids using MPI.          cc
cc                 This routine is always called such that this process       cc
cc                 owns grid "gi1" but does *not* own grid "gi2".             cc
cc                 Also, this routine is always called on some other          cc
cc                 with arguments reversed (obviously) at the same time       cc
cc                 so there is no grid_syncbnd_master.                        cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_syncbnd_slave(gi1, gi2)
      implicit    none
      integer     gi1, gi2
      include     'grid.inc'
      include     'grid_methods.inc'
      include     'glob.inc'
      include     'chr.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'
      integer     parent, nxs,nys,nzs, length
      integer     lx, rx, ly, ry, lz, rz
      integer     nx1, ny1, nz1, nx2, ny2, nz2, i,j,k,
     *            min_i1, max_i1,
     *            min_j1, max_j1,
     *            min_k1, max_k1,
     *            min_i2, max_i2,
     *            min_j2, max_j2,
     *            min_k2, max_k2,
     *            properlength, maxlength,
     *            length1, length2, owner, index
      integer     request_recv, request_send
      logical     first
      logical     fvflag
      real(kind=8) oldnorm1,newnorm1, nrmdiff1, time
      logical     double_equal
      external    double_equal
      integer     calc_bnd_lenb
      external    calc_bnd_lenb
      real(kind=8) myl2norm, myl2norm3d, mynorm
      external     myl2norm, myl2norm3d
      character(3) tmps,tmps2
      integer     gr_chr_gi1
      ! Field number to be output if tracing is on:
      integer     outfnum
      parameter ( outfnum  = 2 )
      logical     blockingsend
      parameter ( blockingsend  = .false. )
      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )
      ! Trace basically just the MPI stuff:
      logical     ltraces
      parameter ( ltraces = .false. )
      ! Just for debugging:
      !integer     s1,s2,s3,s4,s5

      !ltrace  = .false.
      !ltrace2 = .false.
      !if (myid.eq.6) then
      !if (gr_level(gi1).eq.6) then
      !if (gi1 .eq. 161 .or. gi2 .eq. 161) then
         !ltrace  = .true.
         !ltrace2 = .true.
      !end if

 98   format('[',I4,'] ',A,3G15.8)
 99   format('[',I4,'] ',A,4I15)

      if (ltrace) then
         write(*,99) myid,'grid_syncbnd_slave: gi, gj: ', gi1,gi2
         call grid_dump_info(gi1)
         call grid_dump_info(gi2)
      end if

      !
      ! Get pointers to temporary array:
      !
      call load_pointers(gi1)
      gr_chr_gi1 = gr_chr
      owner = grid_return_owner(gi2)

      time = gr_t(gi1)

      nx1 = gr_nx(gi1)
      ny1 = gr_ny(gi1)
      nz1 = gr_nz(gi1)
      call load_scal3d(q(gr_tmp),  0.d0, nx1,ny1,nz1 )
      call load_scal3d(q(gr_flag), 0.d0, nx1,ny1,nz1 )

      nx2 = gr_nx(gi2)
      ny2 = gr_ny(gi2)
      nz2 = gr_nz(gi2)

      call grid_find_intersection(gi1, gi2,length1,
     *            min_i1, max_i1, min_j1, max_j1, min_k1, max_k1)

      !
      ! If no intersection:
      !
      if (length1 .eq. 0) return
      !
      ! If grids do not overlap at their ends
      !         more than the bound_width:
      !
      if (    ((max_i1-min_i1).lt.  bound_width-1 .and.
     *                (max_i1.eq.nx1.or.min_i1.eq.1))
     *    .or.((max_j1-min_j1).lt.  bound_width-1 .and.
     *                (max_j1.eq.ny1.or.min_j1.eq.1))
     *    .or.((max_k1-min_k1).lt.  bound_width-1 .and.
     *                (max_k1.eq.nz1.or.min_k1.eq.1)) ) then
         if (ltrace) then
               write(*,99)myid,'grid_syncbnd_slave:not syncing',gi1,gi2
               call grid_dump_info(gi1)
               call grid_dump_info(gi2)
               write(*,99)myid,'grid_syncbnd_slave:outputting chrgi1'
               call int2str(gi1,tmps)
               call int2str(gi2,tmps2)
             call field_out3d(q(gr_chr),time,     'chrgi1'//tmps//tmps2,
     *             gr_minx(gi1),gr_maxx(gi1),
     *             gr_miny(gi1),gr_maxy(gi1),
     *             gr_minz(gi1),gr_maxz(gi1),
     *             nx1,ny1,nz1,myid)
            write(*,99)myid,'grid_syncbnd_slave: min_i,j,k1= ',min_i1,
     *                      min_j1,min_k1
            write(*,99)myid,'grid_syncbnd_slave: max_i,j,k1= ',max_i1,
     *                      max_j1,max_k1
         end if
         return
      end if

      call grid_find_intersection(gi2, gi1,length2,
     *            min_i2, max_i2, min_j2, max_j2, min_k2, max_k2)


      !
      ! Not needed anymore:
      lx = 1
      rx = 1
      ly = 1
      ry = 1
      lz = 1
      rz = 1

      !
      ! Above lengths don't work for just the boundary:
      !
      length1 = calc_bnd_lenb(min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *                        min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *                        lx,rx, ly,ry, lz,rz,
     *                        nx1,ny1,nz1, nx2,ny2,nz2,bound_width)
      length2 = calc_bnd_lenb(min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *                        min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *                        lx,rx, ly,ry, lz,rz,
     *                        nx2,ny2,nz2, nx1,ny1,nz1,bound_width)


      maxlength = length1
      if (length2 .gt. maxlength) maxlength = length2
      !
      ! If no "effective" intersection:
      !
      if (maxlength .eq. 0) return
      !
      ! If overlap is just one point wide, do not sync
      ! since those points can't be interior
      !
      if (  (min_i1.eq.max_i1 .and .max_i1.eq.nx1) .or.
     *      (min_j1.eq.max_j1 .and .max_j1.eq.ny1) .or.
     *      (min_k1.eq.max_k1 .and .max_k1.eq.nz1) .or.
     *      (min_i2.eq.max_i2 .and .max_i2.eq.nx2) .or.
     *      (min_j2.eq.max_j2 .and .max_j2.eq.ny2) .or.
     *      (min_k2.eq.max_k2 .and .max_k2.eq.nz2) ) then
        if(ltrace)write(*,99)myid,'grid_syncbnd_slv:Notsyncing',gi1,gi2
        if(ltrace)call grid_dump_info(gi1)
        if(ltrace)call grid_dump_info(gi2)
        return
      end if


      if (ltrace.or.(length1.ne.length2)) then
        call int2str(gi1,tmps)
        call int2str(gi2,tmps2)
      write(*,99)myid,'grid_syncbnd_slave: gi1       = ', gi1
      write(*,99)myid,'grid_syncbnd_slave: gi2       = ', gi2
      call field_dump_info(q(gr_chr),nx1,ny1,nz1)
      write(*,99)myid,'grid_syncbnd_slave: length1&2 = ',length1,length2
      write(*,99)myid,'grid_syncbnd_slave: maxlength = ',maxlength
      write(*,99)myid,'grid_syncbnd_slave: owner     = ',owner
      write(*,99)myid,'grid_syncbnd_slave: min_i,j,k1= ',min_i1,min_j1,
     * min_k1
      write(*,99)myid,'grid_syncbnd_slave: max_i,j,k1= ',max_i1,max_j1,
     * max_k1
      write(*,99)myid,'grid_syncbnd_slave: nx,ny,nz1 = ',nx1,ny1,nz1
      write(*,99)myid,'grid_syncbnd_slave: min_i,j,k2= ',min_i2,min_j2,
     * min_k2
      write(*,99)myid,'grid_syncbnd_slave: max_i,j,k2= ',max_i2,max_j2,
     * max_k2
      write(*,99)myid,'grid_syncbnd_slave: nx,ny,nz2 = ',nx2,ny2,nz2
      write(*,99)myid,'grid_syncbnd_slave: lx, rx    = ',lx,rx
      write(*,99)myid,'grid_syncbnd_slave: ly, ry    = ',ly,ry
      write(*,99)myid,'grid_syncbnd_slave: lz, rz    = ',lz,rz
      write(*,98)myid,'grid_syncbnd_slave: gr_t(gi1) = ',time     
       oldnorm1 = myl2norm3d(q(gfunc_pointer(outfnum,gi1)),nx1,ny1,nz1)
       write(*,98)myid,'grid_syncbnd: oldnorm1  = ',oldnorm1
       call field_out3d(q(gfunc_pointer(outfnum,gi1)),time,
     *             'preslave'//tmps//tmps2,
     *             gr_minx(gi1),gr_maxx(gi1),
     *             gr_miny(gi1),gr_maxy(gi1),
     *             gr_minz(gi1),gr_maxz(gi1),
     *             nx1,ny1,nz1,myid)
!     if (gi1.eq.94 .or. gi2.eq.94) then
!       call int2str(gi1,tmps)
!       call int2str(gi2,tmps2)
!       call field_out3d(q(gr_psi4R),gr_t(gi1),
!    *             'pres'//tmps//tmps2,
!    *             gr_minx(gi1),gr_maxx(gi1),
!    *             gr_miny(gi1),gr_maxy(gi1),
!    *             gr_minz(gi1),gr_maxz(gi1),
!    *             nx1,ny1,nz1,myid)
!     end if
c     if (gi1.eq.94 .and. gi2.eq.22) then
      call field_out3d(q(gr_chr),time,      'chrgi1',
     *             gr_minx(gi1),gr_maxx(gi1),
     *             gr_miny(gi1),gr_maxy(gi1),
     *             gr_minz(gi1),gr_maxz(gi1),
     *             nx1,ny1,nz1,myid)
c     else if (gi2.eq.13 .and. gi1.eq.22) then
c     call field_out3d(q(gr_chr),0.d0, 'maskgi2',
c    *             gr_minx(gi1),gr_maxx(gi1),
c    *             gr_miny(gi1),gr_maxy(gi1),
c    *             gr_minz(gi1),gr_maxz(gi1),
c    *             nx1,ny1,nz1,0)
c     end if
      end if

      first = .true.
      do i = 1, num_gfuncs
         !
         ! Keep in mind that this boundary syncing 
         ! takes place *during* the iterative crank-nicholson
         ! and hence we apply the sync'ing to the advanced time level
         !
         if ((gfunc_type(i) .eq.GFUNC_INTEGRAL_ADV)   .or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRAL  .and.(time.eq.0)).or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRAL_DER)
     *                 .or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRAL_DEREL)
     *                 .or.
     *       (gfunc_class(i) .eq. GFUNC_W)
     *                 .or.
     *       (gfunc_name(i) .eq. 'mask')
     *                 .or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRAL_ADVFV).or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRALFV.and.time.eq.0).or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRAL_DERFV) ) then
            !
            if (ltrace2) write(*,99)myid,'  doing field i= ',i
            if (ltraces) write(*,99)myid,'  doing field i= ',i
            !
            if ( gfunc_type(i) .eq. GFUNC_INTEGRALFV .or.
     *           gfunc_type(i) .eq. GFUNC_INTEGRAL_DERFV ) then
               ! Finite Volume fields get sync'ed differently:
                fvflag = .true.
            else
                fvflag = .false.
            end if
            !
            ! Non-blocking recv first:
            !
            if (ltrace2) write(*,99)myid,'  doing IRecv    i= ',i
            if (ltrace2) write(*,98)myid,'  |flag|=',
     *                                     myl2norm(q(gr_flag),length2)
            if (ltraces) write(*,98)myid,'  |flag|=',
     *                                     myl2norm(q(gr_flag),length2)
            call MPI_IRecv(q(gr_flag), maxlength, MPI_DOUBLE_PRECISION,
     *           owner, TAG_SYNC_DATA, MPI_COMM_WORLD,request_recv,ierr)
            if (ltraces) write(*,99)myid,'  IRecv   ierr=',ierr
            if (ltraces) write(*,99)myid,'  IRecv   owner=',owner
            if (ltraces) write(*,99)myid,'  request_recv=',request_recv
            if (ltrace2) write(*,99)myid,'  IRecv   ierr=',ierr
            if (ltrace2) write(*,99)myid,'  request_recv=',request_recv
            !
            ! Make sure tmp space has been sent, so we can write over it:
            !
            if (.not.first) then
               if (.not.blockingsend) then
                  if(ltrace2.or.ltraces)
     *               write(*,99)myid,' Waiting.. send',request_send
                  call MPI_Wait(request_send, status, ierr)
                  if(ltrace2)write(*,99)myid,'  Wait on send done.'
                  if(ltraces)write(*,99)myid,'  Wait on send done.',ierr
               end if
               if(ltrace2)write(*,99)myid,'  Wait ierr: ',ierr
            else
               if(ltraces)write(*,99)myid,'  No wait on first'
               first = .false.
            end if
            !
            ! Copy data into tmp now that it's safe:
            !
            if(ltrace2.or.ltraces)write(*,99)myid,'  Copying data from '
     *                                                   ,gi1
            call copy_int_bndd(   q(gr_tmp),
     *                           q(gfunc_pointer(i,gi1)),
     *                           q(gr_chr_gi1),
     *            min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *            min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *                           lx,rx, ly,ry, lz,rz,
     *                           length1, properlength,
     *            nx1,ny1,nz1, nx2,ny2,nz2, bound_width,fvflag)
            if (properlength.ne.length1) then
               write(*,99)myid,'grid_syncbnd_slave:1ne',length1,
     *                 properlength
               write(*,99)myid,'grid_syncbnd_slave:1  ',length2,
     *                 maxlength
        call int2str(gi1,tmps)
        call int2str(gi2,tmps2)
      write(*,99)myid,'grid_syncbnd_slave: gi1       = ', gi1
      write(*,99)myid,'grid_syncbnd_slave: gi2       = ', gi2
      write(*,99)myid,'grid_syncbnd_slave: properlength',properlength
      write(*,99)myid,'grid_syncbnd_slave: length1&2 = ',length1,length2
      write(*,99)myid,'grid_syncbnd_slave: maxlength = ',maxlength
      write(*,99)myid,'grid_syncbnd_slave: owner     = ',owner
      write(*,99)myid,'grid_syncbnd_slave: min_i,j,k1= ',min_i1,min_j1,
     * min_k1
      write(*,99)myid,'grid_syncbnd_slave: max_i,j,k1= ',max_i1,max_j1,
     * max_k1
      write(*,99)myid,'grid_syncbnd_slave: nx,ny,nz1 = ',nx1,ny1,nz1
      write(*,99)myid,'grid_syncbnd_slave: min_i,j,k2= ',min_i2,min_j2,
     * min_k2
      write(*,99)myid,'grid_syncbnd_slave: max_i,j,k2= ',max_i2,max_j2,
     * max_k2
      write(*,99)myid,'grid_syncbnd_slave: nx,ny,nz2 = ',nx2,ny2,nz2
      write(*,99)myid,'grid_syncbnd_slave: lx, rx    = ',lx,rx
      write(*,99)myid,'grid_syncbnd_slave: ly, ry    = ',ly,ry
      write(*,99)myid,'grid_syncbnd_slave: lz, rz    = ',lz,rz
      write(*,98)myid,'grid_syncbnd_slave: gr_t(gi1) = ',time     
                   write(*,99)myid,'grid_syncbnd_slave: field: ',i
            end if
            if (properlength.ne.length2) then
               write(*,99)myid,'grid_syncbnd_slave:2ne',length2,
     *                 properlength
               write(*,99)myid,'grid_syncbnd_slave:2  ',length1,
     *                 maxlength
        call int2str(gi1,tmps)
        call int2str(gi2,tmps2)
      write(*,99)myid,'grid_syncbnd_slave: gi1       = ', gi1
      write(*,99)myid,'grid_syncbnd_slave: gi2       = ', gi2
      write(*,99)myid,'grid_syncbnd_slave: properlength',properlength
      write(*,99)myid,'grid_syncbnd_slave: length1&2 = ',length1,length2
      write(*,99)myid,'grid_syncbnd_slave: maxlength = ',maxlength
      write(*,99)myid,'grid_syncbnd_slave: owner     = ',owner
      write(*,99)myid,'grid_syncbnd_slave: min_i,j,k1= ',min_i1,min_j1,
     * min_k1
      write(*,99)myid,'grid_syncbnd_slave: max_i,j,k1= ',max_i1,max_j1,
     * max_k1
      write(*,99)myid,'grid_syncbnd_slave: nx,ny,nz1 = ',nx1,ny1,nz1
      write(*,99)myid,'grid_syncbnd_slave: min_i,j,k2= ',min_i2,min_j2,
     * min_k2
      write(*,99)myid,'grid_syncbnd_slave: max_i,j,k2= ',max_i2,max_j2,
     * max_k2
      write(*,99)myid,'grid_syncbnd_slave: nx,ny,nz2 = ',nx2,ny2,nz2
      write(*,99)myid,'grid_syncbnd_slave: lx, rx    = ',lx,rx
      write(*,99)myid,'grid_syncbnd_slave: ly, ry    = ',ly,ry
      write(*,99)myid,'grid_syncbnd_slave: lz, rz    = ',lz,rz
      write(*,98)myid,'grid_syncbnd_slave: gr_t(gi1) = ',time     
                   write(*,99)myid,'grid_syncbnd_slave: field: ',i
            end if
            if (ltrace.and. .false.) then
               !mynorm = myl2norm3d(q(gr_tmp),nx1,ny1,nz1)
               mynorm = myl2norm(q(gr_tmp),properlength)
               if (mynorm.gt.1d7) then
                   write(*,99)myid,'grid_syncbnd_slave: Problem',gi1,gi2
                   write(*,99)myid,'grid_syncbnd_slave: field: ',i
                   write(*,98)myid,'grid_syncbnd_slave: norm: ',mynorm
                   !call field_dump_info(q(gr_tmp),nx1,ny1,nz1)
                 call field_dump_info1D(q(gr_tmp),properlength,'gr_tmp')
               end if
            end if
            !
            ! Non-blocking send:
            !
            if (.not.blockingsend) then
               if(ltraces)write(*,99)myid,'lengths:',length1,length2,
     *                             properlength,maxlength
               if(ltrace2)write(*,99)myid,'  Doing ISend'
               if(ltraces)write(*,98)myid,'  Doing ISend',
     *                       myl2norm(q(gr_tmp),properlength)
             call MPI_ISend(q(gr_tmp),maxlength,MPI_DOUBLE_PRECISION,
     *           owner, TAG_SYNC_DATA,MPI_COMM_WORLD, request_send,ierr)
               if (ltrace2.or.ltraces) then
                   write(*,99)myid,'  ISend ierr: ',ierr
                   write(*,99)myid,'  ISend requestsd:',request_send
               end if
            else
               if (ltrace2) write(*,99)myid,'  Blocking send'
             call MPI_Send(q(gr_tmp),maxlength,MPI_DOUBLE_PRECISION,
     *           owner, TAG_SYNC_DATA,MPI_COMM_WORLD, ierr)
               if (ltrace2) write(*,99)myid,'  MPI_Send ierr=',ierr
            end if
            !
            if (ltrace2) then
               write(*,99)myid,' properlength =',properlength
               write(*,99)myid,' Sending data to  proc: ',owner
               write(*,99)myid,'             of length: ',length1
               write(*,99)myid,'                w/ tag: ',TAG_SYNC_DATA
               write(*,99)myid,'                w/ com: ',MPI_COMM_WORLD
            end if

            !
            ! Wait until data has been received:
            !
            if(ltrace2)write(*,99)myid,' Waiting.. recv',request_recv
            if(ltraces)write(*,99)myid,' Waiting.. recv',request_recv
            call MPI_Wait(request_recv, status, ierr)
            if(ltraces)write(*,99)myid,'  Wait ierr: ',ierr
            if(ltrace2)write(*,99)myid,'  Wait ierr: ',ierr
            if (ltrace2) write(*,99)myid,'  Wait on recv done.'
            !
            if (ltrace2) write(*,99) myid,'  Pasting data from ',owner
            !call load_scal3d(q(gr_flag), -gi2*1.d0, nx1,ny1,nz1 )
            call paste_int_bndd(  q(gr_flag),
     *                           q(gfunc_pointer(i,gi1)),
     *                           q(gr_chr_gi1),
     *            min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *            min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *                           lx,rx, ly,ry, lz,rz,
     *                           nx2,ny2,nz2,
     *                           properlength,
     *                           nx1,ny1,nz1,bound_width,fvflag )
         else
            if(ltrace2)write(*,99)myid,'syncbnd_slave:ignorefield i= ',i
         end if
         !
      end do
      !..................................

      if (.not.first .and. .not.blockingsend) then
         if(ltraces)write(*,99)myid,' Waiting.. send',request_send
         if(ltrace2)write(*,99)myid,' Waiting.. send',request_send
         call MPI_Wait(request_send, status, ierr)
         if(ltrace2)write(*,99)myid,'  Wait ierr: ',ierr
         if(ltraces)write(*,99)myid,'  Wait ierr: ',ierr
      end if

      if (ltrace) then
         write(*,98)myid,'grid_syncbnd_slave: Zeroing gr_flag',
     .                  myl2norm(q(gr_flag),properlength)
         write(*,98)myid,'grid_syncbnd_slave: Zeroing gr_tmp',
     .                  myl2norm(q(gr_tmp),properlength)
      end if
      call load_scal3d(q(gr_flag), 0.d0, nx1,ny1,nz1 )
      call load_scal3d(q(gr_tmp),  0.d0, nx1,ny1,nz1 )


!     if (gi1.eq.94 .or. gi2.eq.94) then
      if (ltrace2) then
       newnorm1 = myl2norm3d(q(gfunc_pointer(outfnum,gi1)),nx1,ny1,nz1)
       nrmdiff1 = abs(newnorm1-oldnorm1)/oldnorm1
       write(*,98)myid,'grid_syncbnd_slave: newnorm1  = ',newnorm1
       write(*,98)myid,'grid_syncbnd_slave: nrmdiff1  = ',nrmdiff1
       if (nrmdiff1 .gt. 0.1) then
          write(*,98)myid,'grid_syncbnd_slave: PROBLEM!!!!! '
          write(*,98)myid,'grid_syncbnd_slave:norms1:',newnorm1,oldnorm1
       end if
       call field_out3d(q(gfunc_pointer(outfnum,gi1)),time,
     *             'postslave'//tmps//tmps2,
     *             gr_minx(gi1),gr_maxx(gi1),
     *             gr_miny(gi1),gr_maxy(gi1),
     *             gr_minz(gi1),gr_maxz(gi1),
     *             nx1,ny1,nz1,myid)
      end if
      if (ltrace.or.ltrace2) then
       call grid_dump_info(gi1)
       !call grid_test(gi1)
         call level_test_local()
       write(*,99) myid,'grid_syncbnd_slave: Done. ',gi1,gi2
      end if
      if(ltrace)then
         write(*,99) myid,'grid_syncbnd_slave: level_test_local',gi1,gi2
         call level_test_local()
         !call load_pointers(gi1)
         write(*,99)myid,'grid_syncbnd_slave: Done:',gi1,gi2
      end if

      return
      end    ! END: grid_syncbnd_slave

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_syncbnd:                                                             cc
cc                 Sync the boundaries of a pair of grids.                    cc
cc                 Where a grid's boundary lies inside the other grid,        cc
cc                 reset those values based on the values of the interior     cc
cc                 the other grid.                                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_syncbnd(gi1, gi2)
      implicit    none
      integer     gi1, gi2
      include     'grid.inc'
      include     'grid_methods.inc'
      include     'glob.inc'
      include     'chr.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'
      integer     parent, nxs,nys,nzs, length1, length2
      integer     lx, rx, ly, ry, lz, rz
      integer     nx1, ny1, nz1, nx2, ny2, nz2, i,j,k,
     *            min_i1, max_i1,
     *            min_j1, max_j1,
     *            min_k1, max_k1,
     *            min_i2, max_i2,
     *            min_j2, max_j2,
     *            min_k2, max_k2,   index,
     *            tmpptr1, tmpptr2, properlength
      logical     fvflag
      logical     double_equal
      external    double_equal
      character(3) tmps,tmps2
      integer     calc_bnd_lenb
      external    calc_bnd_lenb
      real(kind=8) myl2norm3d, mynorm, myl2norm
      real(kind=8) oldnorm1,newnorm1,
     *             oldnorm2,newnorm2,
     *             nrmdiff1,nrmdiff2
      external     myl2norm3d, myl2norm
      integer     gr_chr_gi1, gr_chr_gi2
      ! Field number to be output if tracing is on:
      integer     outfnum
      parameter ( outfnum  = 2 )
      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )

      !if (gr_level(gi1).eq.2 .and. gr_t(gi1).gt.20) then
      !ltrace = .false.
      !ltrace2 = .false.
      !if (gr_level(gi1).eq.1.and.gr_t(gi1).gt.0.0001)then
      !if (gi1 .eq. 161 .or. gi2 .eq.161) then
      !if (myid.eq.6) then
         !ltrace = .true.
         !ltrace2 = .true.
      !end if

 96   format('[',I3,'] ',A,A)
 97   format('[',I3,'] ',A,I5,3F10.5)
 98   format('[',I3,'] ',A,3F10.5)
 99   format('[',I3,'] ',A,3I7)
      !
      ! Get pointers to temporary array:
      !
      call load_pointers(gi2)
      gr_chr_gi2 = gr_chr
      call load_pointers(gi1)
      gr_chr_gi1 = gr_chr


      nx1 = gr_nx(gi1)
      ny1 = gr_ny(gi1)
      nz1 = gr_nz(gi1)

      ! Initialize temporary storage to zero:
      call load_scal3d(q(gr_tmp),  0.d0, nx1,ny1,nz1 )
      call load_scal3d(q(gr_flag), 0.d0, nx1,ny1,nz1 )

      nx2 = gr_nx(gi2)
      ny2 = gr_ny(gi2)
      nz2 = gr_nz(gi2)

      call grid_find_intersection(gi1, gi2, length1,
     *            min_i1, max_i1, min_j1, max_j1, min_k1, max_k1)

      !
      ! If no intersection:
      !
      if (length1 .eq. 0) return
      !
      ! If grids do not overlap at their ends
      !         more than the bound_width:
      !
      if (    ((max_i1-min_i1).lt.bound_width-1 .and.
     *                (max_i1.eq.nx1.or.min_i1.eq.1))
     *    .or.((max_j1-min_j1).lt.bound_width-1 .and.
     *                (max_j1.eq.ny1.or.min_j1.eq.1))
     *    .or.((max_k1-min_k1).lt.  bound_width-1 .and.
     *                (max_k1.eq.nz1.or.min_k1.eq.1)) ) then
          !
         if (ltrace) then
            write(*,99)myid,'grid_syncbnd:not syncing',gi1,gi2
            call grid_dump_info(gi1)
            call grid_dump_info(gi2)
         write(*,99)myid,'grid_syncbnd:min_i,j,k1=',min_i1,min_j1,min_k1
         write(*,99)myid,'grid_syncbnd:max_i,j,k1=',max_i1,max_j1,max_k1
         write(*,99)myid,'grid_syncbnd:nx,ny,nz1 =',nx1,ny1,nz1
         end if
         return
      end if

      call grid_find_intersection(gi2, gi1, length2,
     *            min_i2, max_i2, min_j2, max_j2, min_k2, max_k2)

      !
      ! Not needed anymore:
      lx = 1
      rx = 1
      ly = 1
      ry = 1
      lz = 1
      rz = 1

      !
      ! The above returned values of the length
      ! are for a vector holding the *entire* 
      ! intersection. We only need to store part of the boundary:
      !
      length1 = calc_bnd_lenb(min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *                        min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *                        lx,rx, ly,ry, lz,rz,
     *                        nx1,ny1,nz1, nx2,ny2,nz2,bound_width)
      length2 = calc_bnd_lenb(min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *                        min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *                        lx,rx, ly,ry, lz,rz,
     *                        nx2,ny2,nz2, nx1,ny1,nz1,bound_width)


      if (length1.gt. nx1*ny1*nz1) then
         write(*,99)myid,'grid_syncbnd: Not enough space for g1'
         call my_exit('grid_syncbnd: Not enough space for g1')
      end if
      if (length2.gt. nx2*ny2*nz2) then
         write(*,99)myid,'grid_syncbnd: Not enough space for g2'
         call my_exit('grid_syncbnd: Not enough space for g2')
      end if
      !
      ! If no "effective" intersection:
      !
      if (length1 .eq. 0) return
      !
      ! If overlap is just one point wide, do not sync
      ! since those points can't be interior
      !
      if (  (min_i1.eq.max_i1 .and .max_i1.eq.nx1) .or.
     *      (min_j1.eq.max_j1 .and .max_j1.eq.ny1) .or.
     *      (min_k1.eq.max_k1 .and .max_k1.eq.nz1) .or.
     *      (min_i2.eq.max_i2 .and .max_i2.eq.nx2) .or.
     *      (min_j2.eq.max_j2 .and .max_j2.eq.ny2) .or.
     *      (min_k2.eq.max_k2 .and .max_k2.eq.nz2) ) then
        if(ltrace)write(*,99)myid,'grid_syncbnd: Not syncing ',gi1,gi2
        if(ltrace)call grid_dump_info(gi1)
        if(ltrace)call grid_dump_info(gi2)
        return
      end if

      if (ltrace2) then
         write(*,99) myid,'grid_syncbnd: gi1, gi2  = ', gi1, gi2
         !call grid_dump_info(gi1)
         !call grid_dump_info(gi2)
         write(*,99)myid,'grid_syncbnd: gi1    = ', gi1
         write(*,99)myid,'grid_syncbnd: gi2    = ', gi2
         write(*,99)myid,'grid_syncbnd: min_i1 = ', min_i1
         write(*,99)myid,'grid_syncbnd: max_i1 = ', max_i1
         write(*,99)myid,'grid_syncbnd: min_j1 = ', min_j1
         write(*,99)myid,'grid_syncbnd: max_j1 = ', max_j1
         write(*,99)myid,'grid_syncbnd: min_k1 = ', min_k1
         write(*,99)myid,'grid_syncbnd: max_k1 = ', max_k1
         write(*,99)myid,'grid_syncbnd: min_i2 = ', min_i2
         write(*,99)myid,'grid_syncbnd: max_i2 = ', max_i2
         write(*,99)myid,'grid_syncbnd: min_j2 = ', min_j2
         write(*,99)myid,'grid_syncbnd: max_j2 = ', max_j2
         write(*,99)myid,'grid_syncbnd: min_k2 = ', min_k2
         write(*,99)myid,'grid_syncbnd: max_k2 = ', max_k2
         write(*,99)myid,'grid_syncbnd: nX1    = ', nx1,ny1,nz1
         write(*,99)myid,'grid_syncbnd: nX2    = ', nx2,ny2,nz2
         write(*,99)myid,'grid_syncbnd: length1= ', length1
         write(*,99)myid,'grid_syncbnd: length2= ', length2
      end if
      if (ltrace) then
       !call level_tree_dump()
       write(*,99)myid,'grid_syncbnd: gi1       = ', gi1
       call grid_dump_info(gi1)
       write(*,99)myid,'grid_syncbnd: gi2       = ', gi2
       call grid_dump_info(gi2)
       write(*,99)myid,'grid_syncbnd: length1= ', length1
       write(*,99)myid,'grid_syncbnd: length2= ', length2
       write(*,99)myid,'grid_syncbnd: min_i,j,k1= ',min_i1,min_j1,min_k1
       write(*,99)myid,'grid_syncbnd: max_i,j,k1= ',max_i1,max_j1,max_k1
       write(*,99)myid,'grid_syncbnd: nx,ny,nz1 = ',nx1,ny1,nz1
       write(*,99)myid,'grid_syncbnd: min_i,j,k2= ',min_i2,min_j2,min_k2
       write(*,99)myid,'grid_syncbnd: max_i,j,k2= ',max_i2,max_j2,max_k2
       write(*,99)myid,'grid_syncbnd: nx,ny,nz2 = ',nx2,ny2,nz2
       oldnorm1 = myl2norm3d(q(gfunc_pointer(outfnum,gi1)),nx1,ny1,nz1)
       oldnorm2 = myl2norm3d(q(gfunc_pointer(outfnum,gi2)),nx2,ny2,nz2)
       write(*,98)myid,'grid_syncbnd: oldnorm1  = ',oldnorm1
       write(*,98)myid,'grid_syncbnd: oldnorm2  = ',oldnorm2
       call int2str(gi1,tmps)
       call int2str(gi2,tmps2)
         call field_out3d(q(gfunc_pointer(outfnum,gi1)),gr_t(gi1),
     *             'pres'//tmps//tmps2,
     *             gr_minx(gi1),gr_maxx(gi1),
     *             gr_miny(gi1),gr_maxy(gi1),
     *             gr_minz(gi1),gr_maxz(gi1),
     *             nx1,ny1,nz1,myid)
         call field_out3d(q(gfunc_pointer(outfnum,gi2)),gr_t(gi2),
     *             'pres'//tmps2//tmps,
     *             gr_minx(gi2),gr_maxx(gi2),
     *             gr_miny(gi2),gr_maxy(gi2),
     *             gr_minz(gi2),gr_maxz(gi2),
     *             nx2,ny2,nz2,myid)
      end if

      do i = 1, num_gfuncs
         !
         ! Keep in mind that this boundary syncing 
         ! takes place *during* the iterative crank-nicholson
         ! and hence we apply the sync'ing to the advanced time level
         !
         if ( gfunc_type(i) .eq.GFUNC_INTEGRAL_ADV   .or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRAL  .and.gr_t(gi1).eq.0).or.
     *        gfunc_type(i) .eq.GFUNC_INTEGRAL_DER
     *                 .or.
     *        gfunc_type(i) .eq.GFUNC_INTEGRAL_DEREL
     *                 .or.
     *        gfunc_class(i) .eq. GFUNC_W
     *                 .or.
     *        gfunc_name(i) .eq. 'mask'
     *                 .or.
     *        gfunc_name(i) .eq. 'error'
     *                 .or.
     *        gfunc_type(i) .eq.GFUNC_INTEGRAL_ADVFV .or.
     *       (gfunc_type(i) .eq.GFUNC_INTEGRALFV.and.gr_t(gi1).eq.0).or.
     *        gfunc_type(i) .eq.GFUNC_INTEGRAL_DERFV ) then
            !
            if ( gfunc_type(i) .eq. GFUNC_INTEGRALFV .or.
     *           gfunc_type(i) .eq. GFUNC_INTEGRAL_DERFV ) then
               ! Finite Volume fields get sync'ed differently:
                fvflag = .true.
            else
                fvflag = .false.
            end if
            !
            if (ltrace2)
     *    write(*,96)myid,'grid_syncbnd: gfunc_name:',gfunc_name(i)(1:7)
            !
            call copy_int_bndd(   q(gr_tmp),
     *                           q(gfunc_pointer(i,gi1)),
     *                           q(gr_chr_gi1),
     *            min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *            min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *                           lx,rx, ly,ry, lz,rz,
     *                           length1, properlength,
     *            nx1,ny1,nz1, nx2,ny2,nz2, bound_width,fvflag)
            if (ltrace2) then
               write(*,99)myid,'grid_syncbnd:length1     =',length1
               write(*,99)myid,'grid_syncbnd:properlength=',properlength
               mynorm = myl2norm(q(gr_tmp),properlength)
               write(*,98)myid,'grid_syncbnd:mynorm=',mynorm
               if (mynorm.gt.1d7 .and. .false.) then
                   write(*,99)myid,'grid_syncbnd: Problem',gi1,gi2
                   write(*,98)myid,'grid_syncbnd: norm: ',mynorm
                 call field_dump_info1D(q(gr_tmp),properlength,'gr_tmp')
                 call field_dump_infoB(q(gfunc_pointer(i,gi1)),nx1,ny1,
     .                 nz1,'gr_tmp')
               end if
            end if
            !call load_scal3d(q(gr_tmp), -gi1*1.d0, nx1,ny1,nz1 )
!           Just for debugging:
!           call load_scal_add3d(q(gr_tmp),-gi1*1.d0, 
!    .                     nx1,ny1,nz1)
            call paste_int_bndd(  q(gr_tmp),
     *                           q(gfunc_pointer(i,gi2)),
     *                           q(gr_chr_gi2),
     *            min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *            min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *                           lx,rx, ly,ry, lz,rz,
     *                           nx1,ny1,nz1, properlength,
     *                           nx2,ny2,nz2,bound_width,fvflag )
            !
            !
            call copy_int_bndd(   q(gr_tmp),
     *                           q(gfunc_pointer(i,gi2)),
     *                           q(gr_chr_gi2),
     *            min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *            min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *                           lx,rx, ly,ry, lz,rz,
     *                           length2, properlength,
     *            nx2,ny2,nz2, nx1,ny1,nz1,bound_width,fvflag)
            !call load_scal3d(q(gr_tmp), -gi2*1.d0, nx1,ny1,nz1 )
!           Just for debugging:
!           call load_scal_add3d(q(gr_tmp),-gi2*1.d0,
!    .                      nx1,ny1,nz1)
            if (ltrace2) then
               write(*,99)myid,'grid_syncbnd:length1     =',length1
               write(*,99)myid,'grid_syncbnd:properlength=',properlength
               !mynorm = myl2norm3d(q(gr_tmp),nx1,ny1,nz1)
               mynorm = myl2norm(q(gr_tmp),properlength)
               write(*,98)myid,'grid_syncbnd:mynorm=',mynorm
            end if
            call paste_int_bndd(  q(gr_tmp),
     *                           q(gfunc_pointer(i,gi1)),
     *                           q(gr_chr_gi1),
     *            min_i1,max_i1,min_j1,max_j1,min_k1,max_k1,
     *            min_i2,max_i2,min_j2,max_j2,min_k2,max_k2,
     *                           lx,rx, ly,ry, lz,rz,
     *                           nx2,ny2,nz2, properlength,
     *                           nx1,ny1,nz1,bound_width,fvflag )
            !
            if (ltrace2) then
               write(*,99) myid,'grid_syncbnd: level_test_ afterfield',i
               call level_test_local()
               call load_pointers(gi1)
            end if
            !
         end if
         !
      end do

      ! Large numbers are likely left in the temporary storage and
      ! so we zero them out so that code that monitors for anomolies
      ! does not think there is a real problem:
      call load_scal3d(q(gr_tmp),  0.d0, nx1,ny1,nz1 )
      call load_scal3d(q(gr_flag), 0.d0, nx1,ny1,nz1 )

      if (ltrace .or. ltrace2) then
       newnorm1 = myl2norm3d(q(gfunc_pointer(outfnum,gi1)),nx1,ny1,nz1)
       newnorm2 = myl2norm3d(q(gfunc_pointer(outfnum,gi2)),nx2,ny2,nz2)
       nrmdiff1 = abs(newnorm1-oldnorm1)/oldnorm1
       nrmdiff2 = abs(newnorm2-oldnorm2)/oldnorm2
       write(*,98)myid,'grid_syncbnd: newnorm1  = ',newnorm1
       write(*,98)myid,'grid_syncbnd: newnorm2  = ',newnorm2
       write(*,98)myid,'grid_syncbnd: nrmdiff1  = ',nrmdiff1
       write(*,98)myid,'grid_syncbnd: nrmdiff2  = ',nrmdiff2
       if (nrmdiff1 .gt. 0.1) then 
          write(*,99)myid,'grid_syncbnd: PROBLEM!!!!! '
          write(*,98)myid,'grid_syncbnd: norms1:',newnorm1,oldnorm1
       end if
       if (nrmdiff2 .gt. 0.1) then 
          write(*,98)myid,'grid_syncbnd: PROBLEM!!!!! '
          write(*,98)myid,'grid_syncbnd: norms2:',newnorm2,oldnorm2
       end if
         call field_out3d(q(gfunc_pointer(outfnum,gi1)),gr_t(gi1),
     *             'pos'//tmps//tmps2,
     *             gr_minx(gi1),gr_maxx(gi1),
     *             gr_miny(gi1),gr_maxy(gi1),
     *             gr_minz(gi1),gr_maxz(gi1),
     *             nx1,ny1,nz1,myid)
         call field_out3d(q(gfunc_pointer(outfnum,gi2)),gr_t(gi2),
     *             'pos'//tmps2//tmps,
     *             gr_minx(gi2),gr_maxx(gi2),
     *             gr_miny(gi2),gr_maxy(gi2),
     *             gr_minz(gi2),gr_maxz(gi2),
     *             nx2,ny2,nz2,myid)
         call grid_dump_info(gi1)
         call grid_dump_info(gi2)
         write(*,99)myid,'grid_syncbnd: Calling grid_test.',gi1,gi2
         !call grid_test(gi1)
         !call grid_test(gi2)
         call level_test_local()
         write(*,99)myid,'grid_syncbnd: Done.',gi1,gi2
      end if

      if(ltrace) then
      write(*,99) myid,'grid_syncbnd: level_test_local',gi1,gi2
      call level_test_local()
      end if

      if(ltrace) write(*,99) myid,'grid_syncbnd: Done:',gi1,gi2

      return
      end    ! END: grid_syncbnd


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_inject:                                                              cc
cc                 Update parent grid with values from                        cc
cc                 fine grids. Only injects evolved quantities.               cc
cc                 Allow   derived fields to be computed on parent grid.      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_inject(child,parent)
      implicit    none
      integer     child, parent
      include     'grid.inc'
      include     'glob.inc'
      integer     nxs,nys,nzs, parent_mask_pointer
      integer     nxp, nyp, nzp, nxc, nyc, nzc, i
      integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
      integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
      integer     tmp_ghostwidth
      logical     double_equal
      external    double_equal
      integer     length
      integer     myid, proc_return_myid
      character(3) tmps,tmps2
      ! Field number to be output if tracing is on:
      integer     outfnum
      parameter ( outfnum  = 1 )
      logical     ltrace
      parameter ( ltrace = .false. )

      !ltrace = .false.
      !if (child.eq.179) then
      !ltrace = .true.
      !end if

      !call grid_check_mask(parent,10)

      myid = proc_return_myid()
      !
      ! Store the address of the parent's mask function:
      !
      call load_pointers(parent)
      parent_mask_pointer = gr_chr

      call load_pointers(child)

      nxp = gr_nx(parent)
      nyp = gr_ny(parent)
      nzp = gr_nz(parent)

      nxc = gr_nx(child)
      nyc = gr_ny(child)
      nzc = gr_nz(child)

      !
      ! After an elliptic solve initially, we inject but want
      ! to inject with ghostwidth==1 no matter what the real
      ! value is:
      !
      if (double_equal(gr_t(child),0.d0)) then
         tmp_ghostwidth = 1
      else
         tmp_ghostwidth = ghostwidth
      end if

      !
      ! Find indices which bound intersection region:
      !
      call grid_find_intersection(child, parent, length,
     *            min_ic, max_ic, min_jc, max_jc, min_kc, max_kc)
      call grid_find_intersection(parent, child, length,
     *            min_ip, max_ip, min_jp, max_jp, min_kp, max_kp)

      !
      ! Overlap region on the parent
      !    but with resolution of the parent grid:
      !
      nxs = max_ip - min_ip + 1
      nys = max_jp - min_jp + 1
      nzs = max_kp - min_kp + 1

      if (nxs*nys*nzs .gt. nxc*nyc*nzc) then
         write(*,*) 'grid_inject: Problem. Not enough storage'
         write(*,*) 'grid_inject: nx/y/zs   = ',nxs,nys,nzs
         write(*,*) 'grid_inject: nx/y/zc   = ',nxc,nyc,nzc
         write(*,*) 'grid_inject: nx*y*zs   = ',nxs*nys*nzs
         write(*,*) 'grid_inject: nx*y*zc   = ',nxc*nyc*nzc
         call my_exit('Not enough storage in grid_inject')
      end if

      if (ltrace) then
        call int2str(child,tmps)
        call int2str(parent,tmps2)
         write(*,99)myid,'grid_inject: parent    = ', parent
         write(*,99)myid,'grid_inject: child     = ', child
         write(*,99)myid,'grid_inject: length    = ', length
         write(*,99)myid,'grid_inject: nx/y/zs   = ',nxs,nys,nzs
         write(*,99)myid,'grid_inject: nx/y/zc   = ',nxc,nyc,nzc
         write(*,99)myid,'grid_inject: nx/y/zp   = ',nxp,nyp,nzp
         write(*,99)myid,'grid_inject: min/max_ic= ',min_ic,max_ic
         write(*,99)myid,'grid_inject: min/max_jc= ',min_jc,max_jc
         write(*,99)myid,'grid_inject: min/max_kc= ',min_kc,max_kc
         write(*,99)myid,'grid_inject: min/max_ip= ',min_ip,max_ip
         write(*,99)myid,'grid_inject: min/max_jp= ',min_jp,max_jp
         write(*,99)myid,'grid_inject: min/max_kp= ',min_kp,max_kp
         write(*,98)myid,'grid_inject: time(parent): ',gr_t(parent)
         write(*,98)myid,'grid_inject: time(child) : ',gr_t(child)
         call field_dump_info(q(gfunc_pointer(outfnum,parent)),
     *             nxp,nyp,nzp)
         call field_dump_info(q(gfunc_pointer(outfnum,child)),
     *             nxc,nyc,nzc)
      end if

      call load_scal3d(q(gr_tmp), 0.d0, nxc,nyc,nzc )

      do i = 1, num_gfuncs
         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL .or.
     *        gfunc_name(i) .eq. 'mask'         .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL.or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRAL_DER ) then
            if (i.eq.outfnum .and. ltrace)
     *      call field_out3d(q(gr_chr), gr_t(child),
     *                   'InjectChr'//tmps//tmps2,
     *               gr_minx(child),gr_maxx(child),
     *               gr_miny(child),gr_maxy(child),
     *               gr_minz(child),gr_maxz(child),
     *               nxc,nyc,nzc,myid)
            if (i.eq.outfnum .and. ltrace)
     *      call field_out3d(q(gfunc_pointer(i,child)),gr_t(child),
     *                   'InjectChi'//tmps//tmps2,
     *               gr_minx(child),gr_maxx(child),
     *               gr_miny(child),gr_maxy(child),
     *               gr_minz(child),gr_maxz(child),
     *               nxc,nyc,nzc,myid)
            call restrict_fieldGWb(q(gr_tmp),
     *                           q(gfunc_pointer(i,child)),
     *                           q(gr_chr),
     *                           nxs,nys,nzs,
     *                           nxc,nyc,nzc,
     *                           min_ic,max_ic,
     *                           min_jc,max_jc,
     *                           min_kc,max_kc,
     *                           refine_factor,tmp_ghostwidth)
            !call load_scal1d(q(gr_tmp),-2.d0,nxs*nys*nzs)
            if (i.eq.outfnum .and. ltrace)
     *      call field_out3d(q(gr_tmp),gr_t(parent),
     *           'InjectTmp'//tmps//tmps2,
     *       q(gr_x(parent)+(min_ip-1)),q(gr_x(parent)+(max_ip-1)),
     *       q(gr_y(parent)+(min_jp-1)),q(gr_y(parent)+(max_jp-1)),
     *       q(gr_z(parent)+(min_kp-1)),q(gr_z(parent)+(max_kp-1)),
     *             nxs,nys,nzs,myid)
            call replace_part_field_avgGWc( 
     *                           q(gfunc_pointer(i,parent)),
     *                           q(parent_mask_pointer),
     *                           q(gr_tmp),
     *                           min_ip,
     *                           min_jp,
     *                           min_kp,
     *                           nxp,nyp,nzp,
     *                           nxs,nys,nzs,
     *                           nxc,nyc,nzc,
     *                           min_ic,max_ic,
     *                           min_jc,max_jc,
     *                           min_kc,max_kc,
     *                           refine_factor, tmp_ghostwidth,
     *                           bound_width )
            if (i.eq.outfnum .and. ltrace)
     *      call field_out3d(q(gfunc_pointer(i,parent)),gr_t(parent),
     *                   'InjectPar'//tmps//tmps2,
     *               gr_minx(parent),gr_maxx(parent),
     *               gr_miny(parent),gr_maxy(parent),
     *               gr_minz(parent),gr_maxz(parent),
     *               nxp,nyp,nzp,myid)
         else if ( gfunc_type(i) .eq. GFUNC_INTEGRALFV .or.
     *             gfunc_type(i) .eq. GFUNC_INTEGRAL_DERFV ) then
            call restrict_fieldFV(q(gr_tmp),
     *                           q(gfunc_pointer(i,child)),
     *                           q(gr_chr),
     *                           nxs,nys,nzs,
     *                           nxc,nyc,nzc,
     *                           min_ic,max_ic,
     *                           min_jc,max_jc,
     *                           min_kc,max_kc,
     *                           refine_factor,tmp_ghostwidth)
            !call load_scal1d(q(gr_tmp),-1.d0*child,nxs*nys*nzs)
            call replace_part_field_avgFV(
     *                           q(gfunc_pointer(i,parent)),
     *                           q(parent_mask_pointer),
     *                           q(gr_tmp),
     *                           min_ip,
     *                           min_jp,
     *                           min_kp,
     *                           nxp,nyp,nzp,
     *                           nxs,nys,nzs,
     *                           nxc,nyc,nzc,
     *                           min_ic,max_ic,
     *                           min_jc,max_jc,
     *                           min_kc,max_kc,
     *                           refine_factor, tmp_ghostwidth,
     *                           bound_width )
         !
         end if
      end do

      ! Enforce symmetry conditions after injection:
      if (assume_symmetry.ne.0) call grid_symm(parent)

      if (ltrace)
     *      call field_out3d(q(gfunc_pointer(outfnum,parent)),
     *                   gr_t(parent),
     *                   'InjectPOst'//tmps//tmps2,
     *               gr_minx(parent),gr_maxx(parent),
     *               gr_miny(parent),gr_maxy(parent),
     *               gr_minz(parent),gr_maxz(parent),
     *               nxp,nyp,nzp,myid)

      call inc_mask_overlapB(    q(parent_mask_pointer),
     *                           q(gr_tmp),
     *                           min_ip,
     *                           min_jp,
     *                           min_kp,
     *                           nxp,nyp,nzp,
     *                           nxs,nys,nzs,
     *                           nxc,nyc,nzc,
     *                           min_ic,max_ic,
     *                           min_jc,max_jc,
     *                           min_kc,max_kc,
     *                           refine_factor, tmp_ghostwidth )

      call load_scal3d(q(gr_tmp), 0.d0, nxc,nyc,nzc )

      if (ltrace)
     *     call field_out3d(   q(parent_mask_pointer),
     *             gr_t(parent)*child,'inject_mask',
     *             gr_minx(parent),gr_maxx(parent),
     *             gr_miny(parent),gr_maxy(parent),
     *             gr_minz(parent),gr_maxz(parent),
     *             nxp,nyp,nzp,myid)

      !call grid_check_mask(parent,11)

      if (ltrace) then
         call grid_dump_info(parent)
         call grid_dump_info(child)
         !call grid_test(child)
         call level_test_local()
         write(*,99)myid,'grid_inject: Done.',parent, child
      end if

 97   format('[',I3,'] ',A,I5,3F10.5)
 98   format('[',I3,'] ',A,3F10.5)
 99   format('[',I3,'] ',A,3I5)

      return
      end    ! END: grid_inject

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_intp_bnd_slave:                                                      cc
cc                       Receive boundary data from parent grid               cc
cc                whose process owner runs grid_intp_bnd_master.              cc
cc                Unpacks this data into appropriate grid function            cc
cc                and interpolates in space to provide the Berger &           cc
cc                Oliger boundary update for distributed grids.               cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_intp_bnd_slave(gridnum,parent)
      implicit    none
      integer     gridnum, parent
      include     'grid.inc'
      include     'glob.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'
      integer     nxc, nyc, nzc,
     *            i, length, owner, gwc, nxs, nys, nzs,
     *            bi, ei, bj, ej, bk, ek, nxp, nyp, nzp
      integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
      integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
      real(kind=8)      tp_np1,  tp_n, tc_np1
      real(kind=8)      myl2norm3d, errornorm
      external          myl2norm3d
      logical     grid_return_existence
      character(3) tmps,tmps2
      !
      ! For FV we need our own variables since we do not
      ! want to use an expanded region even if user calls for cubic interp
      !
      integer     min_icFV,max_icFV, min_jcFV,max_jcFV,min_kcFV,max_kcFV
      integer     min_ipFV,max_ipFV, min_jpFV,max_jpFV,min_kpFV,max_kpFV
      integer     nxsFV, nysFV, nzsFV, lengthFV
      integer     c_bnds_C(6), p_bnds_C(6), shp_C(3), length_C, b_C(3)
      integer     c_bnds_L(6), p_bnds_L(6), shp_L(3), length_L, b_L(3)
      integer     c_bnds_Q(6), p_bnds_Q(6), shp_Q(3), length_Q, b_Q(3)


      ! Determine whether it is safe to reset bounds of vfuncs:
      logical     handlevfuncs
      integer     ADV

      ! Field number to be output if tracing is on:
      integer     outfnum
      parameter ( outfnum  = 1 )
      logical     ltrace
      parameter ( ltrace = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(kind=8), allocatable, dimension(:)  :: tdata
      integer :: ntdata, qq
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !ltrace  = .false.
      !ltrace2 = .false.
      !if (gr_level(gridnum).eq.6) then
      !if (gridnum.eq.179) then
      !if (gr_level(gridnum).eq.4.and. myid.eq.3) then
      !if (gr_level(gridnum).eq.4.and. myid.eq.3.and.gridnum.eq.15) then
         !ltrace  = .true.
         !ltrace2 = .true.
      !end if

      owner  = gr_own(parent)

      call grid_find_intersection(gridnum, parent, length,
     *            min_ic, max_ic, min_jc, max_jc, min_kc, max_kc)
      call grid_find_intersection(parent, gridnum, length,
     *            min_ip, max_ip, min_jp, max_jp, min_kp, max_kp)

      ! Keep these for FV fields:
      min_ipFV = min_ip
      min_jpFV = min_jp
      min_kpFV = min_kp
      max_ipFV = max_ip
      max_jpFV = max_jp
      max_kpFV = max_kp
      !
      min_icFV = min_ic
      min_jcFV = min_jc
      min_kcFV = min_kc
      max_icFV = max_ic
      max_jcFV = max_jc
      max_kcFV = max_kc
      !
      nxsFV    = max_ipFV - min_ipFV + 1
      nysFV    = max_jpFV - min_jpFV + 1
      nzsFV    = max_kpFV - min_kpFV + 1
      !
      lengthFV = nxsFV * nysFV * nzsFV

      nxc    = gr_nx(gridnum)
      nyc    = gr_ny(gridnum)
      nzc    = gr_nz(gridnum)

      nxp    = gr_nx(parent)
      nyp    = gr_ny(parent)
      nzp    = gr_nz(parent)

      gwc  = ghostwidth / refine_factor + 1

      !
      ! For higher order methods, in each direction either:
      !     (1) expand the region of the parent you send to the child
      !  or
      !     (2) shrink the region on the child in which you interpolate
      !         ...this shouldn't worry you because in this case some
      !         other parent must overlap in that region anyway.
      !
      !
      !
      if (ltrace) then
      write(*,99)myid,'grid_intp_bnd_slave: Before resetting:'
      write(*,99)myid,'grid_intp_bnd_slavmini/j/kc',min_ic,min_jc,min_kc
      write(*,99)myid,'grid_intp_bnd_slavmaxi/j/kc',max_ic,max_jc,max_kc
      write(*,99)myid,'grid_intp_bnd_slavmini/j/kp',min_ip,min_jp,min_kp
      write(*,99)myid,'grid_intp_bnd_slavmaxi/j/kp',max_ip,max_jp,max_kp
      end if

      call set_interp_bounds(  min_ip, max_ip, min_jp, max_jp,
     &                         min_kp, max_kp, min_ic, max_ic,
     &                         min_jc, max_jc, min_kc, max_kc,
     &                         nxp, nyp, nzp, length,
     &                         bi, bj, bk, refine_factor,
     &                         p_bnds_L, c_bnds_L, shp_L, length_L, b_L,
     &                         p_bnds_C, c_bnds_C, shp_C, length_C, b_C,
     &                         p_bnds_Q, c_bnds_Q, shp_Q, length_Q, b_Q)

      if (c_bnds_Q(2) .lt. c_bnds_Q(1)) return
      if (c_bnds_Q(4) .lt. c_bnds_Q(3)) return
      if (c_bnds_Q(6) .lt. c_bnds_Q(5)) return

      ! People want the quantities computed in analysis
      ! to look "right"...so we can set their boundaries
      ! only if we do so when things are time aligned.
      ! Such quantities do not have two time levels on the
      ! parent so we have to be careful about getting the
      ! data from the right field.
      !
      handlevfuncs = .false.
      if (amrbound_prepost.eq.1.and.amrbound_timealign.eq.0) then
         handlevfuncs = .true.
      end if
      !  Force not to do the vfuncs
      !handlevfuncs = .false.

      if (length_C.gt.nxc*nyc*nzc) then
         write(*,99)myid,'Problem of space in grid_intp_bnd_slave'
      end if

      if (ltrace) then
        call int2str(gridnum,tmps)
        call int2str(parent,tmps2)
         write(*,99)myid,'grid_intp_bnd_slave: Outto: '//tmps//tmps2
      write(*,99)myid,'grid_intp_bnd_slave:     Doing bounds on grid:',
     *                                                   gridnum
      write(*,99)myid,'grid_intp_bnd_slave: gridnum = ', gridnum
      write(*,99)myid,'grid_intp_bnd_slave: parent  = ', parent
      write(*,99)myid,'grid_intp_bnd_slave: owner   = ', owner
      write(*,98)myid,'grid_intp_bnd_slave: time    = ', gr_t(gridnum)
      write(*,99)myid,'grid_intp_bnd_slave: nx/y/zc = ', nxc,nyc,nzc
      write(*,99)myid,'grid_intp_bnd_slave: nx/y/zp = ', nxp,nyp,nzp
      write(*,99)myid,'grid_intp_bnd_slave: gw      = ', ghostwidth
      write(*,99)myid,'grid_intp_bnd_slave: gwc     = ', gwc
      write(*,99)myid,'grid_intp_bnd_slave: length  = ', length
      write(*,99)myid,'grid_intp_bnd_slave: shadow  = ', shadow
      write(*,99)myid,'grid_intp_bnd_slave: nxp*nyp*nzp',nxp*nyp*nzp
      write(*,99)myid,'grid_intp_bnd_slave: amrbound_prepost   = ',
     *                                           amrbound_prepost
      write(*,99)myid,'grid_intp_bnd_slave: amrbound_timealign = ',
     *                                           amrbound_timealign
            if (ltrace2) then
            ADV=0
            call field_out3d(q(gfunc_pointer(outfnum+ADV,gridnum)),
     *             gr_t(gridnum),'prett'//tmps//tmps2,
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             nxc,nyc,nzc,myid)
            ADV=1
            call field_out3d(q(gfunc_pointer(outfnum+ADV,gridnum)),
     *             gr_t(gridnum),'prettADV'//tmps//tmps2,
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             nxc,nyc,nzc,myid)
         end if
      end if

      if (amrbound_prepost.eq.1) call grid_swap_n_np1(gridnum)
      call load_pointers(gridnum)
      call load_scal3d(q(gr_tmp),  0.d0, nxc,nyc,nzc )
      call load_scal3d(q(gr_flag), 0.d0, nxc,nyc,nzc )

      do i = 1, num_gfuncs-1
         ! Set the variables that define the interpolation parameters
         ! for each type of interpolation.
         if (gfunc_interptype(i) .eq. GFUNC_INTERP_LINEAR) then
            call load_interp_const( min_ip, max_ip, min_jp, max_jp,
     &                              min_kp, max_kp, min_ic, max_ic,
     &                              min_jc, max_jc, min_kc, max_kc,
     &                              nxs, nys, nzs, length,
     &                              bi, bj, bk,
     &                              p_bnds_L, c_bnds_L, shp_L,
     &                              length_L, b_L)
         else if (gfunc_interptype(i) .eq. GFUNC_INTERP_CUBIC) then
            call load_interp_const( min_ip, max_ip, min_jp, max_jp,
     &                              min_kp, max_kp, min_ic, max_ic,
     &                              min_jc, max_jc, min_kc, max_kc,
     &                              nxs, nys, nzs, length,
     &                              bi, bj, bk,
     &                              p_bnds_C, c_bnds_C, shp_C,
     &                              length_C, b_C)
         else if (gfunc_interptype(i) .eq. GFUNC_INTERP_QUINTIC) then
            call load_interp_const( min_ip, max_ip, min_jp, max_jp,
     &                              min_kp, max_kp, min_ic, max_ic,
     &                              min_jc, max_jc, min_kc, max_kc,
     &                              nxs, nys, nzs, length,
     &                              bi, bj, bk,
     &                              p_bnds_Q, c_bnds_Q, shp_Q,
     &                              length_Q, b_Q)
         else 
            write(*,99)myid,'grid_intp_bnd_slave: Unknown INTERPTYPE'
            call my_exit('Problem in grid_intp_bnd_slave')
         end if
         !
         if (  gfunc_type(i) .eq. GFUNC_INTEGRAL .or.
     *        (gfunc_type(i) .eq. GFUNC_INTEGRAL_DER .and. handlevfuncs)
     *       ) then
            !
            if (gfunc_type(i).eq.GFUNC_INTEGRAL_DER) then
               ADV = 0
            else
               !
               ! For the u-funcs, put data into the advanced field:
               !
               ADV = 1
            end if
            if(ltrace2) write(*,99)myid,'    ADV = ', ADV
            !
            call MPI_Recv(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_BND_DATA, MPI_COMM_WORLD, status,ierr)
            if(ltrace2) write(*,99)myid,'    Receiving data for ',i
               !    1) Interpolate whole of overlap into work space
               !    2) Paste data from boundaries into field
            if (i.eq.outfnum .and. ltrace)
     *      call field_out3d(q(gr_tmp),gr_t(gridnum),
     *             'tmpsl'//tmps//tmps2,
     *             gr_minx(parent)+(min_ip-1)*gr_h(parent),
     *             gr_minx(parent)+(max_ip-1)*gr_h(parent),
     *             gr_miny(parent)+(min_jp-1)*gr_h(parent),
     *             gr_miny(parent)+(max_jp-1)*gr_h(parent),
     *             gr_minz(parent)+(min_kp-1)*gr_h(parent),
     *             gr_minz(parent)+(max_kp-1)*gr_h(parent),
     *             nxs,nys,nzs,myid)
            if (i.eq.outfnum .and. ltrace) then
      write(*,99)myid,'grid_intp_bnd_slave: nx/y/zs = ', nxs,nys,nzs
      write(*,99)myid,'grid_intp_bnd_slavmini/j/kc',min_ic,min_jc,min_kc
      write(*,99)myid,'grid_intp_bnd_slavmaxi/j/kc',max_ic,max_jc,max_kc
      write(*,99)myid,'grid_intp_bnd_slavmini/j/kp',min_ip,min_jp,min_kp
      write(*,99)myid,'grid_intp_bnd_slavmaxi/j/kp',max_ip,max_jp,max_kp
               write(*,98)myid,'grid_intp_bnd_slave: x:',
     *     q(gr_x(gridnum)+min_ic-1),q(gr_x(gridnum)+max_ic-1)
               write(*,98)myid,'grid_intp_bnd_slave: y:',
     *     q(gr_y(gridnum)+min_jc-1),q(gr_y(gridnum)+max_jc-1)
               write(*,98)myid,'grid_intp_bnd_slave: z:',
     *     q(gr_z(gridnum)+min_kc-1),q(gr_z(gridnum)+max_kc-1)
            end if
            if (gfunc_interptype(i) .eq. GFUNC_INTERP_LINEAR) then
               ! Linear interpolation:
               call init_fieldGWa( q(gr_tmp),
     *                             q(gr_flag),
     *                             nxs,nys,nzs,
     *                             min_ic,min_jc,min_kc,
     *                             max_ic,max_jc,max_kc,
     *                             nxc,nyc,nzc,
     *                             refine_factor,ghostwidth)
            else if (gfunc_interptype(i) .eq. GFUNC_INTERP_CUBIC) then
               ! Cubic interpolation:
#if 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ntdata = nxc*nyc*nzc
               allocate(tdata(ntdata))
               do qq = 1, ntdata
                 tdata(qq) = q(gr_flag + qq -1)
               end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif


               call interp_from_parentB( q(gr_tmp),
     *                          q(gr_flag),
     *                          q(gr_chr), q(gr_mask),
     *                          nxs,nys,nzs,
     *                          min_ic,min_jc,min_kc,
     *                          max_ic,max_jc,max_kc,
     *                          bi,    bj,    bk,
     *                          nxc,nyc,nzc, refine_factor,ghostwidth)
#if 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               call interp_from_parentC( q(gr_tmp),
     *                          tdata,
     *                          q(gr_chr), q(gr_mask),
     *                          nxs,nys,nzs,
     *                          min_ic,min_jc,min_kc,
     *                          max_ic,max_jc,max_kc,
     *                          bi,    bj,    bk,
     *                          nxc,nyc,nzc, refine_factor,ghostwidth)
                call interp_compare(q(gr_flag), tdata,
     *                          min_ic,min_jc,min_kc,
     *                          max_ic,max_jc,max_kc,
     *                          nxc,nyc,nzc)
                deallocate(tdata)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif

            else if (gfunc_interptype(i) .eq. GFUNC_INTERP_QUINTIC) then
               call interp_from_parentC( q(gr_tmp),
     *                          q(gr_flag),
     *                          q(gr_chr), q(gr_mask),
     *                          nxs,nys,nzs,
     *                          min_ic,min_jc,min_kc,
     *                          max_ic,max_jc,max_kc,
     *                          bi,    bj,    bk,
     *                          nxc,nyc,nzc, refine_factor,ghostwidth)
            else if (gfunc_interptype(i) .eq. GFUNC_INTERP_WENO) then
               write(*,*)'WENO interp does not work yet'
            else 
               write(*,*)'Unknown gfunc_interptype'
            end if
!           if (i.eq.1 .and. ltrace)
!    *      call field_out3d(q(gr_flag),gr_t(gridnum),
!    *             'flagpsl'//tmps//tmps2,
!    *             gr_minx(gridnum),gr_maxx(gridnum),
!    *             gr_miny(gridnum),gr_maxy(gridnum),
!    *             gr_minz(gridnum),gr_maxz(gridnum),
!    *             nxc,nyc,nzc,myid)
!           if (i.eq.1 .and. ltrace)
!    *      call field_out3d(q(gr_chr),gr_t(gridnum),
!    *             'chr'//tmps//tmps2,
!    *             gr_minx(gridnum),gr_maxx(gridnum),
!    *             gr_miny(gridnum),gr_maxy(gridnum),
!    *             gr_minz(gridnum),gr_maxz(gridnum),
!    *             nxc,nyc,nzc,myid)
            if (shadow .ne. 0 .and. ADV.eq.1) then
               !
               ! Subtract parent data from data on this grid
               ! to affect the self-shadow hierarchy
               !
               call tre_subtract(     q(gr_flag),
     *                                q(gfunc_pointer(i+1,gridnum)),
     *                                q(gr_error),
     *                                q(gr_chr),
     *                                gfunc_shadow(i),
     *                                min_ic, min_jc, min_kc,
     *                                max_ic, max_jc, max_kc,
     *                                nxc,nyc,nzc, num_injfields,
     *                                refine_factor, ghostwidth)
               if(ltrace)then
                   write(*,97)myid,'grid_intp_bnd_slave: Tre_subtract',
     .                    i,gfunc_shadow(i)
                   write(*,97)myid,'grid_intp_bnd_slave: flag:',i,
     *             myl2norm3d(q(gr_flag),nxc,nyc,nzc)
                   write(*,97)myid,'grid_intp_bnd_slave: ptr :',i,
     *             myl2norm3d(q(gfunc_pointer(i+1,gridnum)),nxc,nyc,nzc)
                   write(*,97)myid,'grid_intp_bnd_slave:com err:',i,
     *             myl2norm3d(q(gr_error),nxc,nyc,nzc)
               end if
            end if
            !
            ! Paste data into boundary regions:
            !
!           if (i.eq.1 .and. ltrace)
!    *      call field_out3d(q(gr_flag),gr_t(gridnum),
!    *             'flagsl'//tmps//tmps2,
!    *             gr_minx(gridnum),gr_maxx(gridnum),
!    *             gr_miny(gridnum),gr_maxy(gridnum),
!    *             gr_minz(gridnum),gr_maxz(gridnum),
!    *             nxc,nyc,nzc,myid)
!           if (i.eq.1 .and. ltrace .and.parent.eq.93)then
!           call paste_boundariesNON(    q(gr_flag),
!    *                                q(gfunc_pointer(i+ADV,gridnum)),
!    *                                q(gr_chr),
!    *                                min_ic, min_jc, min_kc,
!    *                                max_ic, max_jc, max_kc,
!    *                                nxc,nyc,nzc,
!    *                                refine_factor, ghostwidth)
!           end if
            call paste_boundaries(    q(gr_flag),
     *                                q(gfunc_pointer(i+ADV,gridnum)),
     *                                q(gr_chr),
     *                                min_ic, min_jc, min_kc,
     *                                max_ic, max_jc, max_kc,
     *                                nxc,nyc,nzc,
     *                                refine_factor, ghostwidth)
            if (i.eq.outfnum .and. ltrace)
     *      call field_out3d(q(gfunc_pointer(i+ADV,gridnum)),
     *             gr_t(gridnum),'postslADV'//tmps//tmps2,
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             nxc,nyc,nzc,myid)
         else if ( gfunc_type(i) .eq. GFUNC_INTEGRALFV ) then
            !
            ! Finite Volume fields
            !
            ! 1) Interpolate in time on parent       --> gr_tmp
            ! 2) Initialize whole region for FV      --> gr_flag
            ! 3) Paste points into boundary of field --> q()
            call MPI_Recv(q(gr_tmp), lengthFV, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_BND_DATA, MPI_COMM_WORLD, status,ierr)
            call init_fieldFV(    q(gr_tmp),
     *                            q(gr_flag),
     *                            nxsFV,nysFV,nzsFV,
     *                            min_icFV,min_jcFV,min_kcFV,
     *                            max_icFV,max_jcFV,max_kcFV,
     *                            nxc,nyc,nzc, refine_factor,ghostwidth)
            call paste_boundariesFV(  q(gr_flag),
     *                                q(gfunc_pointer(i+1,gridnum)),
     *                                q(gr_chr),
     *                                min_icFV, min_jcFV, min_kcFV,
     *                                max_icFV, max_jcFV, max_kcFV,
     *                                nxc,nyc,nzc,
     *                                refine_factor, ghostwidth)
         end if
      end do
 
      if (shadow.gt.0) then
         if(ltrace)write(*,99)myid,'grid_intp_bnd_slave:Incremnt chr ar'
         call inc_mask_overlapshadow(q(gr_chr),
     *                           q(gr_flag),
     *                           q(gr_error),
     *                           min_ic, min_jc, min_kc,
     *                           max_ic, max_jc, max_kc,
     *                           nxc,nyc,nzc )
            if (ltrace)
     *      call field_out3d(q(gfunc_pointer(outfnum+ADV,gridnum)),
     *             gr_t(gridnum),'postsl2ADV'//tmps//tmps2,
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             nxc,nyc,nzc,myid)
            if (ltrace)
     *      call field_out3d(q(gfunc_pointer(outfnum,gridnum)),
     *             gr_t(gridnum),'postsl2_'//tmps//tmps2,
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             nxc,nyc,nzc,myid)
      end if

      if (amrbound_prepost .eq. 1) then
         call grid_swap_n_np1(gridnum)
      end if

      if (shadow.eq.1) then
         errornorm = myl2norm3d(q(gr_error),nxc,nyc,nzc)
         if (errornorm.gt.1.d8) then
            write(*,99)myid,'grid_intp_bnd_slave: The TRE is HUGE!'
            call field_dump_info(q(gr_error),nxc,nyc,nzc)
      write(*,99)myid,'grid_intp_bnd_slave: gridnum = ', gridnum
      call grid_dump_info(gridnum)
      write(*,99)myid,'grid_intp_bnd_slave: parent  = ', parent
      call grid_dump_info(parent)
      write(*,99)myid,'grid_intp_bnd_slave: owner   = ', owner
      write(*,98)myid,'grid_intp_bnd_slave: time    = ', gr_t(gridnum)
         else if (ltrace) then
            write(*,98)myid,'grid_intp_bnd_slave:Final norm:',errornorm
         end if
      end if

      if (ltrace) then
         write(*,98)myid,'grid_intp_bnd_slave: Zeroing gr_flag:',
     .                     myl2norm3d(q(gr_flag),nxc,nyc,nzc)
         write(*,98)myid,'grid_intp_bnd_slave: Zeroing gr_tmp:',
     .                     myl2norm3d(q(gr_tmp),nxc,nyc,nzc)
      end if
      call load_scal3d(q(gr_flag), 0.d0, nxc,nyc,nzc )
      call load_scal3d(q(gr_tmp),  0.d0, nxc,nyc,nzc )

      if (ltrace) then
         write(*,99)myid,'grid_intp_bnd_slave: grid_test'
         !call grid_test(gridnum)
         call level_test_local()
         write(*,99)myid,'grid_intp_bnd_slave: Done.*****'
      end if

 97   format('[',I3,'] ',A,I5,3F10.5)
 98   format('[',I3,'] ',A,3F10.5)
 99   format('[',I3,'] ',A,3I5)

      return
      end    ! END: grid_intp_bnd_slave

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_2grid_init_master:                                                   cc
cc                         Init a newly created grid using data from a        cc
cc                     nonlocal old sibling grid.                             cc
cc                                                                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_2grid_init_master(grid, sibling)
      implicit    none
      integer     grid, sibling
      include     'grid.inc'
      include     'glob.inc'
      include     'grid_methods.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'
      integer     parent
      integer     nxg, nyg, nzg, nxs, nys, nzs, nxt,nyt,nzt,
     *            i, length, owner,
     *            sowner,
     *            min_i, max_i,
     *            min_j, max_j,
     *            min_k, max_k,
     *            levp, levc, countp, countc
      character(2) tmps,tmps2
      logical     ltrace
      parameter ( ltrace = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )


      owner  = gr_own(grid)

      call load_pointers(sibling)

      nxg    = gr_nx(grid)
      nyg    = gr_ny(grid)
      nzg    = gr_nz(grid)

      nxs    = gr_nx(sibling)
      nys    = gr_ny(sibling)
      nzs    = gr_nz(sibling)

      call grid_find_intersection(sibling,grid, length,
     *            min_i, max_i, min_j, max_j, min_k, max_k)

      if (length.le.0) then
            !
            ! To keep things in sync, need to send zero-length msg
            !
            call MPI_Send(1.d0, 1, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_INIT_DATA, MPI_COMM_WORLD, ierr)
            return
      end if

      nxt = max_i - min_i + 1
      nyt = max_j - min_j + 1
      nzt = max_k - min_k + 1


      if (ltrace) then
 97   format('[',I3,'] ',A,I5,3F10.5)
 98   format('[',I3,'] ',A,3F10.5)
 99   format('[',I3,'] ',A,3I5)
         call int2str(grid,tmps)
         call int2str(sibling,tmps2)
         write(*,99)myid,'grid_2grid_init_master:  grid   = ', grid
         write(*,99)myid,'grid_2grid_init_master: sibling = ', sibling
      write(*,99)myid,'grid_2grid_init_master: level  = ',gr_level(grid)
      write(*,99)myid,'grid_2grid_init_master:levsib=',gr_level(sibling)
         write(*,99)myid,'grid_2grid_init_master:  owner  = ',  owner
         write(*,99)myid,'grid_2grid_init_master: length  = ', length
         write(*,99)myid,'grid_2grid_init_master: nxt     = ', nxt
         write(*,99)myid,'grid_2grid_init_master: nyt     = ', nyt
         write(*,99)myid,'grid_2grid_init_master: nzt     = ', nzt
         write(*,99)myid,'grid_2grid_init_master: nxs     = ', nxs
         write(*,99)myid,'grid_2grid_init_master: nys     = ', nys
         write(*,99)myid,'grid_2grid_init_master: nzs     = ', nzs
         write(*,99)myid,'grid_2grid_init_master: min_i   = ', min_i
         write(*,99)myid,'grid_2grid_init_master: max_i   = ', max_i
         write(*,99)myid,'grid_2grid_init_master: min_j   = ', min_j
         write(*,99)myid,'grid_2grid_init_master: max_j   = ', max_j
         write(*,99)myid,'grid_2grid_init_master: min_k   = ', min_k
         write(*,99)myid,'grid_2grid_init_master: max_k   = ', max_k
           call field_out3d(q(gfunc_pointer(1,sibling)),gr_t(sibling),
     *             'g2gmbefores'//tmps//tmps2,
     *             gr_minx(sibling),gr_maxx(sibling),
     *             gr_miny(sibling),gr_maxy(sibling),
     *             gr_minz(sibling),gr_maxz(sibling),
     *             gr_nx(sibling),gr_ny(sibling),gr_nz(sibling),myid)
      end if

      !
      ! Copies points of the sibling that overlap with grid to
      ! temporary storage
      !
      !
      do i = 1, num_gfuncs
         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL       .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRALFV     .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL .or.
     *        gfunc_name(i) .eq. 'mask'               .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRAL_DERFV .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRAL_DER       ) then
            !
            call copy_overlap(
     *                       q(gr_tmp),
     *                       q(gfunc_pointer(i,sibling)),
     *                       min_i,min_j,min_k,
     *                       nxt,nyt,nzt,
     *                       nxs,nys,nzs )
            if(ltrace2)write(*,99)myid,'grid_2grid_init_master: Sending'
            call MPI_Send(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_INIT_DATA, MPI_COMM_WORLD, ierr)
            if(ltrace2)write(*,99)myid,'grid_2grid_init_master: Sent',i
         end if
      end do
 
      call load_scal3d(q(gr_tmp), 0.d0, nxs,nys,nzs )

      if (ltrace) then
         write(*,99)myid,'grid_2grid_init_master: Done.',grid
      end if

      return
      end    ! END: grid_2grid_init_master

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_2grid_init_slave:                                                    cc
cc                         Init a newly created grid using data from a        cc
cc                     nonlocal old sibling grid.                             cc
cc                                                                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_2grid_init_slave(grid, sibling)
      implicit    none
      integer     grid, sibling
      include     'grid.inc'
      include     'glob.inc'
      include     'grid_methods.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'
      integer     nxg, nyg, nzg, nxt,nyt,nzt,
     *            i, length, sowner,
     *            min_i, max_i,
     *            min_j, max_j,
     *            min_k, max_k,
     *            oi, oj, ok
      character(2) tmps,tmps2
      logical     ltrace
      parameter ( ltrace = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )
      logical     ltracem
      parameter ( ltracem = .false. )


      sowner  = gr_own(sibling)
      !write(*,99) myid, 'grid_2grid_init_slave: sowner: ',sowner

      call load_pointers(grid)

      nxg    = gr_nx(grid)
      nyg    = gr_ny(grid)
      nzg    = gr_nz(grid)

      call grid_find_intersection(sibling,grid, length,
     *            min_i, max_i, min_j, max_j, min_k, max_k)

      if (length.le.0) then
         if (ltracem) write(*,99) myid, '    zero length: ',sowner
            call MPI_Recv(q(gr_tmp), 1, MPI_DOUBLE_PRECISION,
     *               sowner, TAG_INIT_DATA, MPI_COMM_WORLD, status,ierr)
         if (ltracem) write(*,99) myid, '    done', sowner
         return
      end if

      nxt = max_i - min_i + 1
      nyt = max_j - min_j + 1
      nzt = max_k - min_k + 1

      !
      ! Determine offsets between grids:
      !
      if (gr_minx(grid).ge.gr_minx(sibling)) then
         oi = 0
      else
         oi = NINT( (gr_minx(sibling)-gr_minx(grid))/gr_h(grid) )
      end if
      if (gr_miny(grid).ge.gr_miny(sibling)) then
         oj = 0
      else
         oj = NINT( (gr_miny(sibling)-gr_miny(grid))/gr_h(grid) )
      end if
      if (gr_minz(grid).ge.gr_minz(sibling)) then
         ok = 0
      else
         ok = NINT( (gr_minz(sibling)-gr_minz(grid))/gr_h(grid) )
      end if


      if (ltrace2) then
         call int2str(grid,tmps)
         call int2str(sibling,tmps2)
         write(*,99)myid,'grid_2grid_init_slave:  grid   = ', grid
         write(*,99)myid,'grid_2grid_init_slave:  sibling= ', sibling
      end if
      if (ltrace) then
         write(*,99)myid,'grid_2grid_init_slave:  grid   = ', grid
         write(*,99)myid,'grid_2grid_init_slave: sibling = ', sibling
       write(*,99)myid,'grid_2grid_init_slave: level= ',gr_level(grid)
       write(*,99)myid,'grid_2grid_init_slave:levsib=',gr_level(sibling)
         write(*,99)myid,'grid_2grid_init_slave: sowner  = ', sowner
         write(*,99)myid,'grid_2grid_init_slave: length  = ', length
         write(*,99)myid,'grid_2grid_init_slave: nxt     = ', nxt
         write(*,99)myid,'grid_2grid_init_slave: nyt     = ', nyt
         write(*,99)myid,'grid_2grid_init_slave: nzt     = ', nzt
         write(*,99)myid,'grid_2grid_init_slave: oi      = ', oi
         write(*,99)myid,'grid_2grid_init_slave: oj      = ', oj
         write(*,99)myid,'grid_2grid_init_slave: ok      = ', ok
           call field_out3d(q(gfunc_pointer(1,grid)),gr_t(grid),
     *             'g2gsbeforeg'//tmps//tmps2,
     *             gr_minx(grid),gr_maxx(grid),
     *             gr_miny(grid),gr_maxy(grid),
     *             gr_minz(grid),gr_maxz(grid),
     *             gr_nx(grid),gr_ny(grid),gr_nz(grid),myid)
      end if

      !
      ! Copies points of the sibling that overlap with grid to
      ! temporary storage
      !
      !
      if (ltrace2) then
         write(*,99)myid,'grid_2grid_init_slave:Receiving from',sowner
      end if
      if (ltracem) write(*,99) myid, ' Start receiving data...'
      do i = 1, num_gfuncs
         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL       .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRALFV     .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL .or.
     *        gfunc_name(i) .eq. 'mask'               .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRAL_DERFV .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRAL_DER       ) then
            !
            call MPI_Recv(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *               sowner, TAG_INIT_DATA, MPI_COMM_WORLD, status,ierr)
            if (ltracem)
     *         write(*,99) myid, ' recd data ',status(1),status(2),ierr
            call copy_interp_section( 
     *                       q(gfunc_pointer(i,grid)),
     *                       q(gr_tmp),
     *                       oi, oj, ok,
     *                       nxt,nyt,nzt,
     *                       nxg,nyg,nzg, refine_factor)
            if (ltrace2) then
               write(*,99)myid,'grid_2grid_init_slave:       field: ',i
            end if
         end if
      end do
 
      if (ltrace.or.ltrace2) then
           call field_out3d(q(gfunc_pointer(1,grid)),gr_t(grid),
     *             'g2gsafterg'//tmps//tmps2,
     *             gr_minx(grid),gr_maxx(grid),
     *             gr_miny(grid),gr_maxy(grid),
     *             gr_minz(grid),gr_maxz(grid),
     *             gr_nx(grid),gr_ny(grid),gr_nz(grid),myid)
         write(*,99) myid, 'grid_2grid_init_slave: Done.'
      end if

 99   format('[',I3,'] ',A,3I7)

      return
      end    ! END: grid_2grid_init_slave

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_init:                                                                cc
cc                                                                            cc
cc             Initialize a grid from its parent upon creation.               cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_init( child, parent, r)
      implicit    none
      integer     child, parent, r
      include    'glob.inc'
      include    'grid.inc'
      include    'mask.inc'
      integer     nxp, nyp, nzp, nxc, nyc, nzc, i
      integer     nxt, nyt, nzt, length
      integer     bi, bj, bk
      integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
      integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
      integer     parent_mask_ptr
      integer     proc_return_myid, myid
      !
      ! For FV we need our own variables since we do not
      ! want to use an expanded region even if user calls for cubic interp
      !
      integer     min_icFV,max_icFV, min_jcFV,max_jcFV,min_kcFV,max_kcFV
      integer     min_ipFV,max_ipFV, min_jpFV,max_jpFV,min_kpFV,max_kpFV
      integer     nxtFV, nytFV, nztFV, lengthFV
      integer     c_bnds_C(6), p_bnds_C(6), shp_C(3), length_C, b_C(3)
      integer     c_bnds_L(6), p_bnds_L(6), shp_L(3), length_L, b_L(3)
      integer     c_bnds_Q(6), p_bnds_Q(6), shp_Q(3), length_Q, b_Q(3)
      !
      logical     ltrace
      parameter ( ltrace = .false. )
      ! Which field to output if above is true? 
      !    (look in had/include/fields.inc for number)
      integer     fieldnum
      parameter ( fieldnum  = 17 )

      myid = proc_return_myid()

      call load_pointers(parent)
      parent_mask_ptr = gr_mask

      call load_pointers(child)

      nxp = gr_nx(parent)
      nyp = gr_ny(parent)
      nzp = gr_nz(parent)

      nxc = gr_nx(child)
      nyc = gr_ny(child)
      nzc = gr_nz(child)

      call grid_find_intersection(child, parent, length,
     *            min_ic, max_ic, min_jc, max_jc, min_kc, max_kc)
      call grid_find_intersection(parent, child, length,
     *            min_ip, max_ip, min_jp, max_jp, min_kp, max_kp)

      ! Keep these for FV fields:
      min_ipFV = min_ip
      min_jpFV = min_jp
      min_kpFV = min_kp
      max_ipFV = max_ip
      max_jpFV = max_jp
      max_kpFV = max_kp
      !
      min_icFV = min_ic
      min_jcFV = min_jc
      min_kcFV = min_kc
      max_icFV = max_ic
      max_jcFV = max_jc
      max_kcFV = max_kc
      !
      nxtFV    = max_ipFV - min_ipFV + 1
      nytFV    = max_jpFV - min_jpFV + 1
      nztFV    = max_kpFV - min_kpFV + 1
      !
      lengthFV = nxtFV * nytFV * nztFV

      !
      ! For higher order methods, in each direction either:
      !     (1) expand the region of the parent you send to the child
      !  or
      !     (2) shrink the region on the child in which you interpolate
      !         ...this shouldn't worry you because in this case some
      !         other parent must overlap in that region anyway.
      !
      !
      !

      call set_interp_bounds(  min_ip, max_ip, min_jp, max_jp,
     &                         min_kp, max_kp, min_ic, max_ic,
     &                         min_jc, max_jc, min_kc, max_kc,
     &                         nxp, nyp, nzp, length,
     &                         bi, bj, bk, refine_factor,
     &                         p_bnds_L, c_bnds_L, shp_L, length_L, b_L,
     &                         p_bnds_C, c_bnds_C, shp_C, length_C, b_C,
     &                         p_bnds_Q, c_bnds_Q, shp_Q, length_Q, b_Q)
      
      if (c_bnds_Q(2) .lt. c_bnds_Q(1)) return
      if (c_bnds_Q(4) .lt. c_bnds_Q(3)) return
      if (c_bnds_Q(6) .lt. c_bnds_Q(5)) return

      if (length_C.gt.nxc*nyc*nzc) then
         write(*,*) 'Problem of space in grid_init'
         call my_exit('Problem of space in grid_init')
      end if

      if (ltrace) then
 99   format('[',I3,'] ',A,3I7)
         write(*,99)myid,'grid_init: child, parent: ',child,parent
         write(*,99)myid,'grid_init: bi/j/k:        ',b_C(1),b_C(2),
     &                                                b_C(3)
         write(*,99)myid,'grid_init: nx/y/zp:       ',nxp,nyp,nzp
         write(*,99)myid,'grid_init: nx/y/zc:       ',nxc,nyc,nzc
         write(*,99)myid,'grid_init: nx/y/zt:       ',shp_C(1),
     &                                                shp_C(2),shp_C(3)
         write(*,99)myid,'grid_init: min_ic/max_ic: ',c_bnds_C(1),
     &                                                c_bnds_C(2)
         write(*,99)myid,'grid_init: min_jc/max_jc: ',c_bnds_C(3),
     &                                                c_bnds_C(4)
         write(*,99)myid,'grid_init: min_kc/max_kc: ',c_bnds_C(5),
     &                                                c_bnds_C(6)
         write(*,99)myid,'grid_init: level(child):  ',gr_level(child)
         write(*,99)myid,'grid_init: min_ip/max_ip: ',p_bnds_C(1),
     &                                                p_bnds_C(2)
         write(*,99)myid,'grid_init: min_jp/max_jp: ',p_bnds_C(3),
     &                                                p_bnds_C(4)
         write(*,99)myid,'grid_init: min_kp/max_kp: ',p_bnds_C(5),
     &                                                p_bnds_C(6)
         write(*,99)myid,'grid_init: '
      end if

      !
      ! Initialize mask from parent first
      ! since we need it below:
      !
      if (num_masks .gt. 0) then
         ! Use linear interp values for the mask.
         call load_interp_const( min_ip, max_ip, min_jp, max_jp,
     &                            min_kp, max_kp, min_ic, max_ic,
     &                            min_jc, max_jc, min_kc, max_kc,
     &                            nxt, nyt, nzt, length,
     &                            bi, bj, bk,
     &                            p_bnds_L, c_bnds_L, shp_L,
     &                            length_L, b_L)

 
         ! First need to transfer the mask over:
         call copy_overlap( q(gr_tmp),
     *                       q(parent_mask_ptr),
     *                       min_ip, min_jp, min_kp,
     *                       nxt,nyt,nzt,
     *                       nxp,nyp,nzp )
         call paste_mask( q(gr_tmp),
     *                       q(gr_mask),
     *                       nxt,nyt,nzt,
     *                       min_ic,min_jc,min_kc,
     *                       max_ic,max_jc,max_kc,
     *                       bi,    bj,    bk,
     *                       nxc,nyc,nzc,refine_factor)
         if (ltrace) 
     *   call field_out3d(q(gr_mask),gr_t(child),'msk',
     *       gr_minx(child),gr_maxx(child),
     *       gr_miny(child),gr_maxy(child),
     *       gr_minz(child),gr_maxz(child),
     *             nxc,nyc,nzc,proc_return_myid())
      end if

      do i = 1, num_gfuncs
         ! Set the variables that define the interpolation parameters
         ! for each type of interpolation.
         if (gfunc_interptype(i) .eq. GFUNC_INTERP_LINEAR) then
           call load_interp_const(  min_ip, max_ip, min_jp, max_jp,
     &                              min_kp, max_kp, min_ic, max_ic,
     &                              min_jc, max_jc, min_kc, max_kc,
     &                              nxt, nyt, nzt, length,
     &                              bi, bj, bk,
     &                              p_bnds_L, c_bnds_L, shp_L,
     &                              length_L, b_L)
         else if (gfunc_interptype(i) .eq. GFUNC_INTERP_CUBIC) then
           call load_interp_const(  min_ip, max_ip, min_jp, max_jp,
     &                              min_kp, max_kp, min_ic, max_ic,
     &                              min_jc, max_jc, min_kc, max_kc,
     &                              nxt, nyt, nzt, length,
     &                              bi, bj, bk,
     &                              p_bnds_C, c_bnds_C, shp_C,
     &                              length_C, b_C)
         else if (gfunc_interptype(i) .eq. GFUNC_INTERP_QUINTIC) then
           call load_interp_const(  min_ip, max_ip, min_jp, max_jp,
     &                              min_kp, max_kp, min_ic, max_ic,
     &                              min_jc, max_jc, min_kc, max_kc,
     &                              nxt, nyt, nzt, length,
     &                              bi, bj, bk,
     &                              p_bnds_Q, c_bnds_Q, shp_Q,
     &                              length_Q, b_Q)
         end if

         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL       .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRAL_DER ) then
            !
            ! Copy overlap area (perhaps a bit more):
            !
            call copy_overlapb(
     *                       q(gr_tmp),
     *                       q(gfunc_pointer(i,parent)),
     *                       q(parent_mask_ptr),
     *                       min_ip, min_jp, min_kp,
     *                       nxt,nyt,nzt,
     *                       nxp,nyp,nzp )
            if (i.eq.fieldnum .and. ltrace)
     *      call field_out3d(q(gfunc_pointer(i,parent)),gr_t(parent),
     *          'init_par', gr_minx(parent),gr_maxx(parent),
     *                      gr_miny(parent),gr_maxy(parent),
     *                      gr_minz(parent),gr_maxz(parent),
     *             nxp,nyp,nzp,0)
            if (i.eq.fieldnum .and. ltrace)
     *      call field_out3d(q(gr_tmp),gr_t(parent),'init_tmp',
     *       q(gr_x(parent)+(min_ip-1)),q(gr_x(parent)+(max_ip-1)),
     *       q(gr_y(parent)+(min_jp-1)),q(gr_y(parent)+(max_jp-1)),
     *       q(gr_z(parent)+(min_kp-1)),q(gr_z(parent)+(max_kp-1)),
     *             nxt,nyt,nzt,0)
            if (gfunc_interptype(i) .eq. GFUNC_INTERP_LINEAR) then
               !
               ! Linear interpolation:
               !
               call init_fieldGWa( q(gr_tmp),
     *                          q(gfunc_pointer(i,child)),
     *                          nxt,nyt,nzt,
     *                          min_ic,min_jc,min_kc,
     *                          max_ic,max_jc,max_kc,
     *                          nxc,nyc,nzc, refine_factor,ghostwidth)
            else if (gfunc_interptype(i) .eq. GFUNC_INTERP_CUBIC) then
               !
               ! Cubic interpolation:
               !
               call interp_from_parentB( q(gr_tmp),
     *                          q(gfunc_pointer(i,child)),
     *                          q(gr_chr), q(gr_mask),
     *                          nxt,nyt,nzt,
     *                          min_ic,min_jc,min_kc,
     *                          max_ic,max_jc,max_kc,
     *                          bi,    bj,    bk,
     *                          nxc,nyc,nzc, refine_factor,ghostwidth)
            else if (gfunc_interptype(i) .eq. GFUNC_INTERP_QUINTIC) then
               !
               ! 5th order interpolation:
               !
               call interp_from_parentC( q(gr_tmp),
     *                          q(gfunc_pointer(i,child)),
     *                          q(gr_chr), q(gr_mask),
     *                          nxt,nyt,nzt,
     *                          min_ic,min_jc,min_kc,
     *                          max_ic,max_jc,max_kc,
     *                          bi,    bj,    bk,
     *                          nxc,nyc,nzc, refine_factor,ghostwidth)
            else
              write(0,*)'Unkonwn interpolation type.'
            end if
            if (i.eq.fieldnum .and. ltrace)
     *      call field_out3d(q(gfunc_pointer(i,child)),gr_t(child),
     *        'init_chld',gr_minx(child),gr_maxx(child),
     *                    gr_miny(child),gr_maxy(child),
     *                    gr_minz(child),gr_maxz(child),
     *                    nxc,nyc,nzc,0)
         else if ( gfunc_type(i) .eq. GFUNC_INTEGRALFV .or.
     *             gfunc_type(i) .eq. GFUNC_INTEGRAL_DERFV ) then
            !
            ! Copy overlap area:
            !
            call copy_overlapb(
     *                       q(gr_tmp),
     *                       q(gfunc_pointer(i,parent)),
     *                       q(parent_mask_ptr),
     *                       min_ipFV, min_jpFV, min_kpFV,
     *                       nxtFV,nytFV,nztFV,
     *                       nxp,nyp,nzp )
            call init_fieldFV(  q(gr_tmp),
     *                          q(gfunc_pointer(i,child)),
     *                          nxtFV,nytFV,nztFV,
     *                          min_icFV,min_jcFV,min_kcFV,
     *                          max_icFV,max_jcFV,max_kcFV,
     *                          nxc,nyc,nzc, refine_factor,ghostwidth)
            !
         end if
      end do

      if (ltrace) then
         write(*,99)myid,'grid_init: Done w/ child/parent:',child,parent
      end if

      return
      end    ! END: grid_init


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_init_parent:                                                         cc
cc                         Copy grid function data within the region          cc
cc                     covered by child to tmp space and send to child.       cc
cc                     Input "grid" is the child grid to be initialized.      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_init_parent(grid,parent)
      implicit    none
      integer     grid, parent
      include     'grid.inc'
      include     'glob.inc'
      include     'mask.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'
      integer     nxp, nyp, nzp, nxc, nyc, nzc, nxt,nyt,nzt,
     *            ax,  ay,  az,  i, length, owner,
     *            levp, levc
      integer     bi, bj, bk
      integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
      integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
      !
      ! For FV we need our own variables since we do not
      ! want to use an expanded region even if user calls for cubic interp
      !
      integer     min_icFV,max_icFV, min_jcFV,max_jcFV,min_kcFV,max_kcFV
      integer     min_ipFV,max_ipFV, min_jpFV,max_jpFV,min_kpFV,max_kpFV
      integer     nxtFV, nytFV, nztFV, lengthFV
      integer     c_bnds_C(6), p_bnds_C(6), shp_C(3), length_C, b_C(3)
      integer     c_bnds_L(6), p_bnds_L(6), shp_L(3), length_L, b_L(3)
      integer     c_bnds_Q(6), p_bnds_Q(6), shp_Q(3), length_Q, b_Q(3)

      logical     ltrace
      parameter ( ltrace = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )


      owner  = gr_own(grid)

      call load_pointers(parent)

      nxp    = gr_nx(parent)
      nyp    = gr_ny(parent)
      nzp    = gr_nz(parent)

      nxc    = gr_nx(grid)
      nyc    = gr_ny(grid)
      nzc    = gr_nz(grid)

      call grid_find_intersection(grid, parent, length,
     *            min_ic, max_ic, min_jc, max_jc, min_kc, max_kc)
      call grid_find_intersection(parent, grid, length,
     *            min_ip, max_ip, min_jp, max_jp, min_kp, max_kp)

      ! Keep these for FV fields:
      min_ipFV = min_ip
      min_jpFV = min_jp
      min_kpFV = min_kp
      max_ipFV = max_ip
      max_jpFV = max_jp
      max_kpFV = max_kp
      !
      min_icFV = min_ic
      min_jcFV = min_jc
      min_kcFV = min_kc
      max_icFV = max_ic
      max_jcFV = max_jc
      max_kcFV = max_kc
      !
      nxtFV    = max_ipFV - min_ipFV + 1
      nytFV    = max_jpFV - min_jpFV + 1
      nztFV    = max_kpFV - min_kpFV + 1
      !
      lengthFV = nxtFV * nytFV * nztFV

      ! For higher order methods, in each direction either:
      !     (1) expand the region of the parent you send to the child
      !  or
      !     (2) shrink the region on the child in which you interpolate
      !         ...this shouldn't worry you because in this case some
      !         other parent must overlap in that region anyway.
      !
      !
      !

      call set_interp_bounds(  min_ip, max_ip, min_jp, max_jp,
     &                         min_kp, max_kp, min_ic, max_ic,
     &                         min_jc, max_jc, min_kc, max_kc,
     &                         nxp, nyp, nzp, length,
     &                         bi, bj, bk, refine_factor,
     &                         p_bnds_L, c_bnds_L, shp_L, length_L, b_L,
     &                         p_bnds_C, c_bnds_C, shp_C, length_C, b_C,
     &                         p_bnds_Q, c_bnds_Q, shp_Q, length_Q, b_Q)

      if (c_bnds_Q(2) .lt. c_bnds_Q(1)) return
      if (c_bnds_Q(4) .lt. c_bnds_Q(3)) return
      if (c_bnds_Q(6) .lt. c_bnds_Q(5)) return

      !
      ! Dimensions of part to be copied:
      !
      if (length_C.gt.nxp*nyp*nzp) then
         write(*,*) 'Problem of space in grid_init_parent'
         call my_exit('Problem of space in grid_init_parent')
      end if

      if (ltrace) then
 99   format('[',I3,'] ',A,3I7)
         write(*,99)myid,'grid_init_parent:  grid   = ', grid
         write(*,99)myid,'grid_init_parent:  owner  = ', owner
         write(*,99)myid,'grid_init_parent: parent  = ', parent
         write(*,99)myid,'grid_init_parent: powner  = ', gr_own(parent)
         write(*,99)myid,'grid_init_parent: length  = ', length
         write(*,99)myid,'grid_init_parent: bi/j/k:        ',bi,bj,bk
         write(*,99)myid,'grid_init_parent: nx/y/zp:       ',nxp,nyp,nzp
         write(*,99)myid,'grid_init_parent: nx/y/zc:       ',nxc,nyc,nzc
         write(*,99)myid,'grid_init_parent: nx/y/zt:       ',nxt,nyt,nzt
         write(*,99)myid,'grid_init_parent: '
         !call grid_dump_info(grid)
         !call grid_dump_info(parent)
      end if

      ! Initialize mask from parent first
      ! since we need it below:
      !
      if (num_masks .gt. 0) then
         ! First need to transfer the mask over:
         call load_interp_const( min_ip, max_ip, min_jp, max_jp,
     &                              min_kp, max_kp, min_ic, max_ic,
     &                              min_jc, max_jc, min_kc, max_kc,
     &                              nxt, nyt, nzt, length,
     &                              bi, bj, bk,
     &                              p_bnds_L, c_bnds_L, shp_L,
     &                              length_L, b_L)

         call copy_overlap( q(gr_tmp),
     *                       q(gr_mask),
     *                       min_ip, min_jp, min_kp,
     *                       nxt,nyt,nzt,
     *                       nxp,nyp,nzp )
         if(ltrace2)write(*,99)myid,'grid_init_parent: Call MPI_Send'
         call MPI_Send(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_INIT_DATA, MPI_COMM_WORLD, ierr)
         if(ltrace2)write(*,99)myid,'grid_init_parent: MPI_Send',owner
      end if

      !
      ! Copies points of the parent that overlap with child to
      ! temporary storage
      !
      !
      do i = 1, num_gfuncs
         if (gfunc_interptype(i) .eq. GFUNC_INTERP_LINEAR) then
           call load_interp_const( min_ip, max_ip, min_jp, max_jp,
     &                              min_kp, max_kp, min_ic, max_ic,
     &                              min_jc, max_jc, min_kc, max_kc,
     &                              nxt, nyt, nzt, length,
     &                              bi, bj, bk,
     &                              p_bnds_L, c_bnds_L, shp_L,
     &                              length_L, b_L)
         else if (gfunc_interptype(i) .eq. GFUNC_INTERP_CUBIC) then
           call load_interp_const( min_ip, max_ip, min_jp, max_jp,
     &                              min_kp, max_kp, min_ic, max_ic,
     &                              min_jc, max_jc, min_kc, max_kc,
     &                              nxt, nyt, nzt, length,
     &                              bi, bj, bk,
     &                              p_bnds_C, c_bnds_C, shp_C,
     &                              length_C, b_C)
         else if (gfunc_interptype(i) .eq. GFUNC_INTERP_QUINTIC) then
           call load_interp_const( min_ip, max_ip, min_jp, max_jp,
     &                              min_kp, max_kp, min_ic, max_ic,
     &                              min_jc, max_jc, min_kc, max_kc,
     &                              nxt, nyt, nzt, length,
     &                              bi, bj, bk,
     &                              p_bnds_Q, c_bnds_Q, shp_Q,
     &                              length_Q, b_Q)
         end if

         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRAL_DER ) then
            if(ltrace2)call
     *         field_dump_info(q(gfunc_pointer(i,parent)),nxp,nyp,nzp)
            call copy_overlapb(
     *                       q(gr_tmp),
     *                       q(gfunc_pointer(i,parent)),
     *                       q(gr_mask),
     *                       min_ip, min_jp, min_kp,
     *                       nxt,nyt,nzt,
     *                       nxp,nyp,nzp )
            if(ltrace2)write(*,99)myid,'grid_init_parent: sending',i
            call MPI_Send(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_INIT_DATA, MPI_COMM_WORLD, ierr)
         else if ( gfunc_type(i) .eq. GFUNC_INTEGRALFV .or.
     *             gfunc_type(i) .eq. GFUNC_INTEGRAL_DERFV ) then
            !
            ! Copy overlap area:
            !
            call copy_overlapb(
     *                       q(gr_tmp),
     *                       q(gfunc_pointer(i,parent)),
     *                       q(gr_mask),
     *                       min_ipFV, min_jpFV, min_kpFV,
     *                       nxtFV,nytFV,nztFV,
     *                       nxp,nyp,nzp )
            if(ltrace2)call
     *         field_dump_info(q(gfunc_pointer(i,parent)),nxp,nyp,nzp)
            if(ltrace2)write(*,99)myid,'grid_init_parentFVsending',i
            call MPI_Send(q(gr_tmp), lengthFV, MPI_DOUBLE_PRECISION,
     *                 owner, TAG_INIT_DATA, MPI_COMM_WORLD, ierr)
         end if
      end do
 
      if (ltrace) then
         write(*,99)myid,'grid_init_parent: Done.'
      end if

      return
      end    ! END: grid_init_parent

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_init_child:                                                          cc
cc                         Copy grid funcation data received from parent      cc
cc                     to grid and interpolate information in order to        cc
cc                     initialize a grid.                                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_init_child(grid, parent)
      implicit    none
      integer     grid, parent
      include     'grid.inc'
      include     'glob.inc'
      include     'mask.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'
      integer     nxc, nyc, nzc, nxt,nyt,nzt, nxp,nyp,nzp,
     *            i, length, powner,
     *            levp, levc, countp, countc, grid_return_level
      integer     bi, bj, bk
      integer     min_ic,max_ic, min_jc,max_jc, min_kc,max_kc
      integer     min_ip,max_ip, min_jp,max_jp, min_kp,max_kp
      logical     grid_return_existence
      !
      ! For FV we need our own variables since we do not
      ! want to use an expanded region even if user calls for cubic interp
      !
      integer     min_icFV,max_icFV, min_jcFV,max_jcFV,min_kcFV,max_kcFV
      integer     min_ipFV,max_ipFV, min_jpFV,max_jpFV,min_kpFV,max_kpFV
      integer     nxtFV, nytFV, nztFV, lengthFV
      integer     c_bnds_C(6), p_bnds_C(6), shp_C(3), length_C, b_C(3)
      integer     c_bnds_L(6), p_bnds_L(6), shp_L(3), length_L, b_L(3)
      integer     c_bnds_Q(6), p_bnds_Q(6), shp_Q(3), length_Q, b_Q(3)

      logical     ltrace
      parameter ( ltrace = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )


      powner = gr_own(parent)

      call load_pointers(grid)

      nxc    = gr_nx(grid)
      nyc    = gr_ny(grid)
      nzc    = gr_nz(grid)

      nxp    = gr_nx(parent)
      nyp    = gr_ny(parent)
      nzp    = gr_nz(parent)

      call grid_find_intersection(grid, parent, length,
     *            min_ic, max_ic, min_jc, max_jc, min_kc, max_kc)
      call grid_find_intersection(parent, grid, length,
     *            min_ip, max_ip, min_jp, max_jp, min_kp, max_kp)

      ! Keep these for FV fields:
      min_ipFV = min_ip
      min_jpFV = min_jp
      min_kpFV = min_kp
      max_ipFV = max_ip
      max_jpFV = max_jp
      max_kpFV = max_kp
      !
      min_icFV = min_ic
      min_jcFV = min_jc
      min_kcFV = min_kc
      max_icFV = max_ic
      max_jcFV = max_jc
      max_kcFV = max_kc
      !
      nxtFV    = max_ipFV - min_ipFV + 1
      nytFV    = max_jpFV - min_jpFV + 1
      nztFV    = max_kpFV - min_kpFV + 1
      !
      lengthFV = nxtFV * nytFV * nztFV

      call set_interp_bounds(  min_ip, max_ip, min_jp, max_jp,
     &                         min_kp, max_kp, min_ic, max_ic,
     &                         min_jc, max_jc, min_kc, max_kc,
     &                         nxp, nyp, nzp, length,
     &                         bi, bj, bk, refine_factor,
     &                         p_bnds_L, c_bnds_L, shp_L, length_L, b_L,
     &                         p_bnds_C, c_bnds_C, shp_C, length_C, b_C,
     &                         p_bnds_Q, c_bnds_Q, shp_Q, length_Q, b_Q)

      if (c_bnds_Q(2) .lt. c_bnds_Q(1)) return
      if (c_bnds_Q(4) .lt. c_bnds_Q(3)) return
      if (c_bnds_Q(6) .lt. c_bnds_Q(5)) return

      if (ltrace) then
 99   format('[',I3,'] ',A,3I7)
         write(*,99)myid,'grid_init_child:  grid   = ', grid
         write(*,99)myid,'grid_init_child: parent  = ', parent
         write(*,99)myid,'grid_init_child: powner  = ', powner
         write(*,99)myid,'grid_init_child: length  = ', length_C
         write(*,99)myid,'grid_init_child: bi/j/k:        ',bi,bj,bk
         write(*,99)myid,'grid_init_child: nx/y/zp:       ',nxp,nyp,nzp
         write(*,99)myid,'grid_init_child: nx/y/zc:       ',nxc,nyc,nzc
         write(*,99)myid,'grid_init_child: nx/y/zt:       ',nxt,nyt,nzt
         write(*,99)myid,'grid_init_child: '
         !call grid_dump_info(grid)
         !call grid_dump_info(parent)
      end if

      !
      ! Initialize mask from parent first
      ! since we need it below:
      !
      if (num_masks .gt. 0) then
         call load_interp_const( min_ip, max_ip, min_jp, max_jp,
     &                            min_kp, max_kp, min_ic, max_ic,
     &                            min_jc, max_jc, min_kc, max_kc,
     &                            nxt, nyt, nzt, length,
     &                            bi, bj, bk,
     &                            p_bnds_L, c_bnds_L, shp_L,
     &                            length_L, b_L)

         if(ltrace2)write(*,99)myid,'grid_init_child:Call mpi_recv'
         call MPI_Recv(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *               powner, TAG_INIT_DATA, MPI_COMM_WORLD, status,ierr)
         call paste_mask( q(gr_tmp),
     *                       q(gr_mask),
     *                       nxt,nyt,nzt,
     *                       min_ic,min_jc,min_kc,
     *                       max_ic,max_jc,max_kc,
     *                       bi,    bj,    bk,
     *                       nxc,nyc,nzc,refine_factor)
         if (ltrace)
     *   call field_out3d(q(gr_mask),gr_t(grid),'msk',
     *       gr_minx(grid),gr_maxx(grid),
     *       gr_miny(grid),gr_maxy(grid),
     *       gr_minz(grid),gr_maxz(grid),
     *             nxc,nyc,nzc,myid)
         if(ltrace2)write(*,99)myid,'grid_init_child:recv mask:',powner
      end if

      !
      ! Copies points of the parent that overlap with child to
      ! temporary storage
      !
      !
      do i = 1, num_gfuncs
         ! Set the variables that define the interpolation parameters
         ! for each type of interpolation.
         if (gfunc_interptype(i) .eq. GFUNC_INTERP_LINEAR) then
           call load_interp_const( min_ip, max_ip, min_jp, max_jp,
     &                              min_kp, max_kp, min_ic, max_ic,
     &                              min_jc, max_jc, min_kc, max_kc,
     &                              nxt, nyt, nzt, length,
     &                              bi, bj, bk,
     &                              p_bnds_L, c_bnds_L, shp_L,
     &                              length_L, b_L)
         else if (gfunc_interptype(i) .eq. GFUNC_INTERP_CUBIC) then
           call load_interp_const( min_ip, max_ip, min_jp, max_jp,
     &                              min_kp, max_kp, min_ic, max_ic,
     &                              min_jc, max_jc, min_kc, max_kc,
     &                              nxt, nyt, nzt, length,
     &                              bi, bj, bk,
     &                              p_bnds_C, c_bnds_C, shp_C,
     &                              length_C, b_C)
         else if (gfunc_interptype(i) .eq. GFUNC_INTERP_QUINTIC) then
           call load_interp_const( min_ip, max_ip, min_jp, max_jp,
     &                              min_kp, max_kp, min_ic, max_ic,
     &                              min_jc, max_jc, min_kc, max_kc,
     &                              nxt, nyt, nzt, length,
     &                              bi, bj, bk,
     &                              p_bnds_Q, c_bnds_Q, shp_Q,
     &                              length_Q, b_Q)
         end if


         if ( gfunc_type(i) .eq. GFUNC_INTEGRAL .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRAL_DEREL .or.
     *        gfunc_type(i) .eq. GFUNC_INTEGRAL_DER ) then
            !
       if(ltrace2)write(*,99)myid,'grid_init_child: recv from:',powner,i
            call MPI_Recv(q(gr_tmp), length, MPI_DOUBLE_PRECISION,
     *               powner, TAG_INIT_DATA, MPI_COMM_WORLD, status,ierr)
            !
            if (gfunc_interptype(i).eq.GFUNC_INTERP_LINEAR) then
               !
               ! Linear interpolation:
               !
               call init_fieldGWa( q(gr_tmp),
     *                          q(gfunc_pointer(i,grid)),
     *                          nxt,nyt,nzt,
     *                          min_ic,min_jc,min_kc,
     *                          max_ic,max_jc,max_kc,
     *                          nxc,nyc,nzc, refine_factor,ghostwidth)
            else if (gfunc_interptype(i).eq.GFUNC_INTERP_CUBIC) then
               !
               ! Cubic interpolation:
               !
               call interp_from_parentb( q(gr_tmp),
     *                          q(gfunc_pointer(i,grid)),
     *                          q(gr_chr), q(gr_mask),
     *                          nxt,nyt,nzt,
     *                          min_ic,min_jc,min_kc,
     *                          max_ic,max_jc,max_kc,
     *                          bi,    bj,    bk,
     *                          nxc,nyc,nzc, refine_factor,ghostwidth)
            else if (gfunc_interptype(i).eq.GFUNC_INTERP_QUINTIC) then
               !
               ! 5th order interpolation:
               !
               call interp_from_parentC( q(gr_tmp),
     *                          q(gfunc_pointer(i,grid)),
     *                          q(gr_chr), q(gr_mask),
     *                          nxt,nyt,nzt,
     *                          min_ic,min_jc,min_kc,
     *                          max_ic,max_jc,max_kc,
     *                          bi,    bj,    bk,
     *                          nxc,nyc,nzc, refine_factor,ghostwidth)
 
            else
               write(*,*)'Unknown gfunc_interptype'
            end if
         else if ( gfunc_type(i) .eq. GFUNC_INTEGRALFV .or.
     *             gfunc_type(i) .eq. GFUNC_INTEGRAL_DERFV ) then
       if(ltrace2)write(*,99)myid,'grid_init_child:FVrecv from',powner,i
            call MPI_Recv(q(gr_tmp), lengthFV, MPI_DOUBLE_PRECISION,
     *               powner, TAG_INIT_DATA, MPI_COMM_WORLD, status,ierr)
            call init_fieldFV(  q(gr_tmp),
     *                          q(gfunc_pointer(i,grid)),
     *                          nxtFV,nytFV,nztFV,
     *                          min_icFV,min_jcFV,min_kcFV,
     *                          max_icFV,max_jcFV,max_kcFV,
     *                          nxc,nyc,nzc, refine_factor,ghostwidth)
         end if
      end do
 
      if (ltrace) then
         !call grid_dump_info(grid)
         write(*,99)myid,'grid_init_child: Done w/grid: ',grid
      end if

      return
      end    ! END: grid_init_child

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  hyperamr_return_level_from_res:                                           cc
cc                  Returns level of grid based upon resolution               cc
cc              NB: this routine is for use on the hyper side of things.      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function hyperamr_return_level_from_res(mydx)
      implicit    none
      real*8      mydx
      include    'param.inc'
      logical     ltrace
      parameter ( ltrace = .true. )

      if(ltrace)write(*,*)'hyperamr_return_level_from_res: mydx=',
     .                    mydx,hxyz0
      if (mydx.le.0) then
         hyperamr_return_level_from_res = -1
      else
         hyperamr_return_level_from_res = log(hxyz0/mydx)/log(2.d0)
      end if
      if (ltrace) write(*,*) 'hyperamr_return_level_from_res: level =',
     .                             hyperamr_return_level_from_res

      return
      end       ! END: hyperamr_return_level_from_res

