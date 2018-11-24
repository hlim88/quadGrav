cc
cc   These routines apply to levels and do *not* include
cc   the user-defined include files grid.inc, fields.inc or mem.inc
cc
cc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_tree_test:                                                          cc
cc            Rudimentary test  of tree structure.                            cc
cc            Only called for debugging purposes.                             cc
cc            NB: Now checks that all grids are ALIVE.                        cc
cc            NB: Now checks that the parents of all grids have a non-null    cc
cc                           child pointer.                                   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_tree_test()
      implicit none
      integer  gi,     lev_i
      integer  parent, child, sibling,
     *         level_return_start,
     *         grid_return_parent,
     *         grid_return_sibling,
     *         grid_return_child,
     *         grid_return_level,
     *         proc_return_myid
      logical  grid_return_existence,
     *         level_grid_exists_within,
     *         grid_alive

      logical     ltrace
      parameter ( ltrace = .false. )

      lev_i = 0
      if (ltrace) write(*,*) '                         Tree Test:'
  2   gi = level_return_start(lev_i)
      if (grid_return_existence(gi) ) then
  3      continue
         parent  = grid_return_parent(gi)
         child   = grid_return_child(gi)
         sibling = grid_return_sibling(gi)
         if (ltrace) write(*,*) ' Checking grid: ',gi
         !
         ! Check ALIVE status:
         !
         if (.not.grid_alive(gi)) then
            write(*,*) ' Grid ',gi,' is not alive'
            call my_exit('Dead grid')
         end if
         !
         ! Check parent:
         !
         if (.not.grid_return_existence(parent).and.(lev_i.gt.0)) then
            write(*,*) '   Grid ',gi,' has nonexistent parent',parent
            call my_exit('Bad parent')
         else if (lev_i.gt.0) then
            if ( .not.level_grid_exists_within(parent, lev_i-1) ) then
               call level_tree_dump()
               write(*,*) '   Grid ',gi,' has parent ',parent,
     *                    ' not in correct level', lev_i-1
               write(*,*) '   myid        = ',proc_return_myid()
               write(*,*) '   parent grid = ', parent
               write(*,*) '   return level= ', grid_return_level(parent)
               write(*,*) '   parent in level?',
     *        level_grid_exists_within(parent,grid_return_level(parent))
               call my_exit('Parent not in proper level')

            else if (.not.
     *          grid_return_existence(grid_return_child(parent))) then
               write(*,*) ' Parent of grid claims not to have children'
               call my_exit('Childless parent in tree_test')
            end if
         end if
         !
         ! Check child:
         !
         if (grid_return_existence(child)) then
            if ( .not. level_grid_exists_within(child, lev_i+1) ) then
               write(*,*) '  Grid ',gi,' has child ',child,
     *                    ' not in correct level'
               call my_exit('Child not in proper level')
            end if
         end if
         !
         ! Check sibling:
         !
         if (grid_return_existence(sibling)) then
            if ( .not. level_grid_exists_within(sibling, lev_i) ) then
               write(*,*) '  Grid ',gi,' has sibling not in correct lev'
               call my_exit('sibling not in proper level')
            end if
         end if
         !
         ! continue with other grids
         !
         gi = grid_return_sibling(gi)
         if (grid_return_existence(gi))  goto 3
         lev_i = lev_i + 1
         goto 2
      end if
      if (ltrace) write(*,*) '                    END: Tree Test:'

      !
      ! See what grid numbers are free:
      !
c     call grid_dump_freenums()

      return
      end      ! END: subroutine level_tree_test

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_time_aligned:                                                       cc
cc            Compares two levels to see if they are possibly                 cc
cc            time aligned for use with outputting grids.                     cc
cc            Presumes that the levels are close to being time aligned        cc
cc            so that it needs only look at the level counters to see         cc
cc            if they are definitely not time aligned.                        cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function level_time_aligned(level1, level2)
      implicit none
      integer  level1, level2
      include 'glob.inc'
      integer  l1, l2

      if (level1 .gt. level2) then
         l1 = level1
         l2 = level2
      else if (level2 .gt. level1) then
         l1 = level2
         l2 = level1
      else 
         !
         ! If a level is compared to itself, definitely time aligned:
         !
         level_time_aligned = .true.
         return
      end if

  10  level_time_aligned = mod(lev_count(l1),refine_factor).eq.0
      l1 = l1 - 1
      if (level_time_aligned .and. l1.ne.l2) goto 10

      return
      end           ! END: level_time_aligned

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_grid_exists_within                                                  cc
cc            Returns true if given grid can be found                         cc
cc            by traversing tree structure.                                   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function level_grid_exists_within(gi,level)
      implicit none
      integer  gi, level
      include 'glob.inc'
      integer  tmp, grid_return_sibling
      logical  grid_return_existence

      level_grid_exists_within = .false.

      tmp = Levelp(level)
   5  if ( tmp .eq. gi ) then
         level_grid_exists_within = .true.
         return
      else if (grid_return_existence(tmp)) then
         tmp = grid_return_sibling(tmp)
         goto 5
      else 
         return
      end if

      return
      end           ! END: level_grid_exists_within

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_dump_count:                                                         cc
cc            Dump the counters for the various levels.                       cc
cc            Diagnostic called from my_exit() to show                        cc
cc            which levels were brought into existence                        cc
cc            during the evolution.                                           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_dump_count()
      implicit none
      include 'mpif.h'
      include 'mpi_stuff.inc'
      include 'glob.inc'
      integer  i

      do i = 0, allowedl
         if (lev_count(i).eq.0) return
         write(*,10) myid, i, lev_count(i)
      end do

 10   format('[',I3,'] level_dump_count:   Level: ',I3,'  Count:',I7)

      return
      end           ! END: level_dump_count

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_return_existence:                                                   cc
cc                  Returns true if the level:                                cc
cc                     --- is between 0 and maxlev                            cc
cc                     --- has at least one grid on that level                cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function level_return_existence(lev)
      implicit    none
      integer     lev
      include    'glob.inc'
      integer     level_return_start   
      logical     grid_return_existence
      external    grid_return_existence

      level_return_existence = .false.

      !
      ! First need to check the level less than
      ! maxlev (Remember than there's a level 0):
      !
      if (lev .lt. maxlev .and. lev .ge. 0) then
         !
         ! Then check to see if there is a grid
         ! at this level:
         !
         if (grid_return_existence(level_return_start(lev)))
     *      level_return_existence = .true.
      end if

      return
      end       ! END: level_return_existence

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_find_first_child:                                                   cc
cc                  Returns the grid number of the first grid to have         cc
cc                  the input grid as its parent.                             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function level_find_first_child(gi)
      implicit    none
      integer     gi
      include    'glob.inc'
      integer     tmpgi
      logical     grid_return_existence
      integer     level_return_start,
     *            grid_return_level,
     *            grid_return_parent,
     *            grid_return_sibling

      !
      ! Set to null grid:
      !
      level_find_first_child = -1

      tmpgi = level_return_start(grid_return_level(gi)+1)

  10  if ( grid_return_existence(tmpgi) ) then
         if ( grid_return_parent(tmpgi) .eq. gi ) then
            level_find_first_child = tmpgi
         else
            tmpgi = grid_return_sibling(tmpgi)
            goto 10
         end if
      end if


      return
      end       ! END: level_find_first_child

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_return_finest_ever:                                                 cc
cc                  Returns the level of finest level that ever existed.      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function level_return_finest_ever()
      implicit    none
      include    'glob.inc'

      level_return_finest_ever = 0

 11   if ( level_return_finest_ever .eq. maxlev-1 ) return

c     write(*,*) 'level_return_finest_ever:', level_return_finest_ever

      if ( lev_count(level_return_finest_ever+1) .gt. 0 ) then
         level_return_finest_ever = level_return_finest_ever + 1
         go to 11
      end if

      return
      end       ! END: level_return_finest_ever

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_return_finest:                                                      cc
cc                  Returns the level of finest existing level.               cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function level_return_finest( Levelp, maxlev )
      implicit    none
      integer     maxlev, Levelp(0:maxlev)

      level_return_finest = 0

 11   if ( level_return_finest .eq. maxlev-1 ) return

      if ( Levelp(level_return_finest+1) .gt. 0 ) then
         level_return_finest = level_return_finest + 1
         go to 11
      end if

      return
      end       ! END: level_return_finest

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_set_mask:                                                           cc
cc            Sets the mask functions (chr()) of grids this processor         cc
cc            owns on this level such that boundaries of these grids          cc
cc            which occur *within* other grids at this level are marked       cc
cc            as CHR_deco_bdy (see chr.inc).                                  cc
cc     NB:  Loops over all siblings, even those with different parents, now.  cc
cc     NB:  Assumes grids can only overlap with other grids with the same     cc
cc          parent (in other words, grids completely contained w/in their     cc
cc          unique parents). This could be relaxed later if necessary.        cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_set_mask(level)
      implicit none
      integer  level
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      gi, parent, lastgrid, gj, lstparent

      logical      ltrace
      parameter (  ltrace = .false. )

      if (ltrace) then
         write(*,99)myid,'level_set_mask: Setting mask on level:',level
      end if

      !
      ! Loop first over all grids and reset their chr() arrays
      !    (we use grid_mk_alive() to re-init chr() since they're all alive anyway)
      !
      gi = level_return_start(level)
  11  if ( grid_return_existence(gi) ) then
         call grid_mk_alive(gi)
         gi = grid_return_sibling(gi)
         goto 11
      end if

      !
      ! Loop over all grids on level:
      !
      gi = level_return_start(level)
  10  if ( grid_return_existence(gi) ) then
         if (ltrace) write(*,99)myid,'     Work w/ grid: ',gi
         gj = grid_return_sibling(gi)
  15     if (grid_return_existence(gj) ) then
            if ( grid_is_local(gi) ) then
               if (ltrace) write(*,99)myid,'        set_mask:',gi,gj
               call grid_set_mask(gi, gj)
            end if
            if (grid_is_local(gj)) then
               if (ltrace) write(*,99)myid,'        set_mask:',gj,gi
               call grid_set_mask(gj, gi)
            end if
            gj = grid_return_sibling(gj)
            goto 15
         end if
         gi = grid_return_sibling(gi)
         goto 10
      end if

      if (ltrace) then
         write(*,99) myid, 'level_set_mask: Deco boundaries set'
         write(*,99) myid, 'level_set_mask: Now cleanup AMR bdy'
      end if

      !
      ! Loop over all grids on level:
      !
      gi = level_return_start(level)
  20  if ( grid_return_existence(gi) ) then
         if (ltrace) write(*,99)myid,'     Work w/ grid: ',gi
         if ( grid_is_local(gi) ) then
            if (ltrace) write(*,99)myid,'        set_amrmask:',gi
            call grid_set_amrmask(gi)
         end if
         gi = grid_return_sibling(gi)
         goto 20
      end if

      if (ltrace) then
         write(*,99) myid, 'level_set_mask: Done'
      end if

 99   format('[',I3,'] ',A,2I5)
      return
      end           ! END: level_set_mask

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_tree_dump:                                                          cc
cc            Dump the tree structure working from the                        cc
cc            pointers to the first grids for each level, Levelp().           cc
cc            Purely for diagnostic/statistical purposes. Does not            cc
cc            change any pointers or other data structures.                   cc
cc                                                                            cc
cc            Modified for MPI: Added output of grid owner.                   cc
cc                              Outputs which processor is writing.           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_tree_dump()
      implicit none
      include 'grid_methods.inc'
      integer  proc_return_myid
      external proc_return_myid
      integer  gi,     num,   lev_i, myid, nx,ny,nz, numspaces
      real(kind=8)      levwork,levh

      integer        indent
      parameter    ( indent = 3 )
      character(256)  spaces
      data           spaces / ' ' /

      logical        outputbrief
      parameter    ( outputbrief= .true. )

      logical        outputdims
      parameter    ( outputdims = .false. )

      logical        ltrace
      parameter    ( ltrace = .false. )


      myid  = proc_return_myid()
      lev_i = 0
      num   = 0
      write(*,*) '                         Tree Dump:'
      if (ltrace) write(*,*) 'Consider      level: ',lev_i
  2   if (level_return_existence(lev_i)) then
         gi        = level_return_start(lev_i)
         numspaces = 1 + lev_i * indent
         if (ltrace) write(*,*) 'Level starts with gi:',gi
         if (outputdims .and. .not.outputbrief) then
            call grid_get_dims(gi,nx,ny,nz)
            write(*,11)  myid,
     *                lev_i,         
     *                spaces(1:numspaces),       gi,
     *                grid_return_sibling(gi),   grid_return_owner(gi),
     *                grid_return_work(gi),      nx,
     *                ny,                        nz,
     *                grid_return_minx(gi),      grid_return_miny(gi),
     *                grid_return_minz(gi),
     *                grid_return_time(gi)
         else if (.not.outputbrief) then
            if (ltrace) write(*,*) 'Outputting...'
            write(*,10)  myid,
     *                lev_i,        
     *                spaces(1:numspaces),       gi,
     *                grid_return_sibling(gi),   grid_return_owner(gi),
     *                grid_return_work(gi),      grid_return_h(gi)
            if (ltrace) write(*,*) 'Done...'
         else
            levwork  = grid_return_work(gi)
            levh     = grid_return_h(gi)
         end if
  5      gi = grid_return_sibling(gi)
         if (grid_return_existence(gi)) then
            num = num + 1
            if (outputdims .and. .not.outputbrief) then
               call grid_get_dims(gi,nx,ny,nz)
               write(*,21) myid,
     *                spaces(1:numspaces),       gi,
     *                grid_return_sibling(gi),   grid_return_owner(gi),
     *                grid_return_work(gi),      nx,
     *                ny,                        nz,
     *                grid_return_minx(gi),      grid_return_miny(gi),
     *                grid_return_minz(gi),
     *                grid_return_time(gi)
            else if (.not.outputbrief) then
               write(*,20) myid,
     *                spaces(1:numspaces),       gi,
     *                grid_return_sibling(gi),   grid_return_owner(gi),
     *                grid_return_work(gi)
            else
               levwork  = levwork  + grid_return_work(gi)
            end if
            goto 5
         end if
         if (outputbrief) then
            write(*,54) '       Level: ',lev_i,'  num grids: ',num+1,
     *              '  work: ',levwork, ' h: ',levh
 54         format(A,I4,A,I4,A,F9.2,A,F10.3)
         end if
         lev_i = lev_i + 1
         if (ltrace) write(*,*) 'Consider next level: ',lev_i
         if (ltrace) write(*,*) 'Consider next level: ',lev_i
         num   = 0
         goto 2
      end if
      write(*,*) '                    END: Tree Dump:'
      call mem_stat_brief()
      !call proc_dump_wkld()


 10   format('[',I4,'] Level ',I2, ' --->',A,I3,'(S:',
     *       I3,' O:',I4,' W:',F7.1,')', ' h:',F5.1)
 20   format('[',I4,']              ',A,I3,'(S:'
     *       ,I3,' O:',I4,' W:',F7.1,')')

 11   format('[',I4,'] Level ',I2, ' --->',A,I3,'(S:',
     *      I3,' O:',I4,' W:',F7.1,')','D(',I3,',',I3,',',I3,') ',
     *      'A(',F7.2,',',F7.2,',',F7.2,')T',F8.4)
 21   format('[',I4,']              ',A,I3,'(S:'
     *      ,I3,' O:',I4,' W:',F7.1,')','D(',I3,',',I3,',',I3,') ',
     *      'A(',F7.2,',',F7.2,',',F7.2,')T',F8.4)

 30   format(A,3F9.3)
 31   format(A,3I5)
      return
      end      ! END: subroutine level_tree_dump

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_return_start                                                        cc
cc            Returns the first grid in the given level.                      cc
cc            It is not guaranteed to exist though.                           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function level_return_start(level)
      implicit none
      integer  level
      include 'glob.inc'

      if (level .ge. maxlev .or. level .lt. 0) then
         level_return_start = -1
      else
         level_return_start = Levelp(level)
      end if

      return
      end           ! END: level_return_start


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_return_work:                                                        cc
cc             Compute "work" (relative to coarse grid) to                    cc
cc             evolve the current level (including all sub-levels).           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(kind=8) function level_return_work( level )
      implicit none
      integer  level
      integer  lev_i,  gr_i
      real(kind=8)      grid_return_work
      integer     grid_return_sibling,
     *            level_return_start
      logical     grid_return_existence,
     *            level_return_existence
      external    grid_return_work,
     *            grid_return_sibling,
     *            level_return_start,
     *            grid_return_existence,
     *            level_return_existence
      logical     ltrace
      parameter ( ltrace = .false. )
     

      if ( .not. level_return_existence(level)   ) then
         write(*,*) 'level_return_work: Level nonexistent, level=',level
         level_return_work = -1.d0
         return
      end if

      if (ltrace) then
         write(*,*) 'level_return_work:      level = ', level
      end if
 
      level_return_work  = 0.d0
      lev_i              = level

      !
      ! Loop over level:
      !
  10  gr_i       = level_return_start(lev_i)
  11  if (grid_return_existence(gr_i)) then
         level_return_work = level_return_work + grid_return_work(gr_i)
         if (ltrace) write(*,*) ' grid, work ',gr_i,level_return_work
         gr_i = grid_return_sibling(gr_i)
         if (ltrace) write(*,*) ' trying sibling   ',gr_i
         goto 11
      end if
      !
      ! Go to next level:
      !
      lev_i      = lev_i + 1
      if (ltrace) write(*,*) ' trying new level ',lev_i
      if ( .not.level_return_existence(lev_i) ) goto 12
      goto 10
 
  12  continue

      if (ltrace) write(*,*) ' Exitting level_return_work'

      return
      end     ! END: function level_return_work

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_inc_count:                                                          cc
cc            Increments the step counter for a given level.                  cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_inc_count(level)
      implicit none
      integer  level
      include 'glob.inc'
      include 'grid_methods.inc'
      integer  gi, myid
      integer  proc_return_myid
      external proc_return_myid
      real(kind=8) time,parenttime

      logical     ltrace
      parameter ( ltrace   = .false. )

      logical     ltrace2
      parameter ( ltrace2  = .false. )

      myid = proc_return_myid()

      if (level .lt. maxlev .and. level .ge. 0) then
         lev_count(level) = lev_count(level) + 1
      else
         write(*,*)'level_inc_count: Problem level does not exist',level
         call my_exit('Cannot increment counter of nonexistent level')
      end if

      if (ltrace2) write(*,99) myid,'level_inc_count: ',
     *       ' Increment counter for level: ',level,lev_count(level)

      if (level .eq. 0 .and.  flush_period.gt.0) then
         if (mod( lev_count(level),flush_period).eq.0) then
            ! every coarse grid step, flush filesystem buffers
            ! on all procs
            !    Does not seem to do the trick on ranger:
            !    and caused a crash on QB for Jan:
            !call system('sync')
            !    Try closing the SDF files (just make sure
            !     that GFT_APPEND is set)
            call gft_close_all()
            if (myid.eq.0) then
               write(*,*)'level_inc_count: Flushing SDF buffers'
               write(*,*)'Environment variable GFT_APPEND should be set'
            end if
         end if
      end if

      !
      ! Loop over all non-local grids and advance their time:
      !
      gi = level_return_start(level)
  10  if ( grid_return_existence(gi) ) then
         if ( .not.grid_is_local(gi) ) then
            if (ltrace)write(*,98)myid,'level_inc_count:',
     *       ' Increment time on grid: ',gi,grid_return_time(gi)
            call grid_advance_time(gi)
         end if
         gi = grid_return_sibling(gi)
         goto 10
      end if

      !
      ! Miguel & Matt have pointed out a problem
      ! in which, at late times, two levels get output to
      ! SDF with differences in times of order 10e-14 and
      ! DV therefore treats them at different times.
      ! This appears to be "simply" an output issue related
      ! to the fact that the grid time is computed by
      ! adding Delta T for each time step. On finer levels,
      ! truncation error will presumably accumulate faster
      ! because of more steps.
      !
      if (level.gt.0) then
         if (level_time_aligned(level,level-1) ) then
            parenttime = grid_return_time(level_return_start(level-1))
            if(ltrace)write(*,97)myid,'level_inc_count',level,parenttime
            gi = level_return_start(level)
  20        if ( grid_return_existence(gi) ) then
               time = grid_return_time(gi)
               call grid_set_time(gi,parenttime)
               if(ltrace )then
                  write(*,97)myid,'level:',gi,time,grid_return_time(gi),
     *                    time-grid_return_time(gi)
               end if
               gi = grid_return_sibling(gi)
               goto 20
            end if
         end if
      end if

 99   format('[',I3,'] ',A,A,2I7)
 98   format('[',I3,'] ',A,A,I7,2F10.4)
 97   format('[',I3,'] ',A,I7,2F10.4,F20.14)

      return
      end           ! END: level_inc_count

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_apply:                                                              cc
cc               Driver to apply various actions                              cc
cc               to all grids in a level.                                     cc
cc               Loops through Levelp()-->gr_sibling() pointers.              cc
cc               Actions listed in file: action.inc                           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_apply(lev,action)
      implicit    none
      integer     lev, action
      include    'action.inc'
      include    'grid_methods.inc'
      include    'mpif.h'
      include    'mpi_stuff.inc'
      include    'glob.inc'
      include    'param.inc'
      include    'tracers.inc'
      integer     gi, owner, powner, sowner, nextproc
      integer     gj, own, level, parent,gnx,gny,gnz,ax,ay,az
      integer     from, iters, j, flev, tmpp
      real(kind=8)starttime, endtime

      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )
      logical     ltraces
      parameter ( ltraces = .false. )


      starttime = MPI_Wtime()
      if (ltrace)write(*,*) 'level_apply: lev, action: ',lev,action

      !
      ! grab first grid on level:
      !
      gi    = level_return_start(lev)
      if ( .not.grid_return_existence(gi) .and.
     *     action.ne.REFINE .and. action.ne.READ_STATE) then
         if (ltrace) write(*,*) 'level_apply: If level does not exist'
         if (ltrace) write(*,*) 'level_apply: then nothing to do'
         return
      end if

      flev = level_return_finest(Levelp, maxlev)

      !.....................................
      !
      !  These actions apply to a level so execute here.
      !  Other actions apply to all grids on a level.
      !
      !.....................................
      if (action .eq. LEVINC) then
         !
         ! Increment level counter:
         !
         if (ltrace)write(*,*) 'level_apply: SendAll LEVINC  ',lev
         call send_actionAll(LEVINC,gi,0,lev,0,0,0,0,0,0,0)
         call level_inc_count(lev)
         return
      else if (action .eq. STEP) then
         !.....................................
         !
         !  With overlapping grids, taking a step involves
         !  the entire level:                              
         !          1) Iterate once on all grids on this level
         !          2) For all boundary regions located in interior of other
         !             grids at same level, replace with those inter. values
         !          3) Repeat (1-2) until converged.
         !          4) Do any necessary smoothing, switch n and np1 levels, etc
         !
         !.....................................
         if(ltraces)write(*,*) 'level_apply: Begin STEP: ',lev
         if (ltrace)write(*,*) 'level_apply: SendAll STEP    ',lev
         call send_actionAll(STEP,gi,0,lev,0,0,0,0,0,0,0)
         do iters = 1, num_evol_iters
            if(ltrace )write(*,*) 'level_apply: STEP:  iters= ',iters
            if(ltraces)write(*,*) 'level_apply: STEP:  iters= ',iters
            call level_iter_local(lev)
            if(ltrace )write(*,*) 'level_apply: STEP:  syncbnd'
            if(ltraces)write(*,*) 'level_apply: STEP:  syncbnd'
            call level_syncbnd_local(lev)
         end do
         if (lev .eq. level_return_finest(Levelp,maxlev) .and.
     *       use_mask .gt. 0  ) then
            if (ltrace) write(*,*) 'level_apply: Calling bh_mask:'
            if (ltraces)write(*,*) 'level_apply: Calling bh_mask:'
            call level_bhmask_local(lev)
            call level_syncbnd_local(lev)
         end if
         if(ltrace )write(*,*) 'level_apply: STEP:   STEP2'
         if(ltraces)write(*,*) 'level_apply: STEP:   STEP2'
         call level_step2_local(lev)
         if(ltrace )write(*,*) 'level_apply: STEP:   Done '
         if(ltraces)write(*,*) 'level_apply: STEP:   Done '
         return
      else if (action .eq. SYNCBND) then
         if (ltrace)write(*,*) 'level_apply: SendAll SYNCBND ',lev
         call send_actionAll(SYNCBND,gi,0,lev,0,0,0,0,0,0,0)
         call level_syncbnd_local(lev)
         return
      else if (action .eq. SYNCBNDI) then
         if (ltrace)write(*,*) 'level_apply: SendAll SYNCBNDI',lev
         call send_actionAll(SYNCBNDI,gi,0,lev,0,0,0,0,0,0,0)
         call level_syncbndi_local(lev)
         return
      else if (action .eq. COMPUT) then
         if (ltrace)write(*,*) 'level_apply: SendAll COMPUT  ',lev
         call send_actionAll(COMPUT,gi,0,lev,0,0,0,0,0,0,0)
         call level_comput_local(lev)
         !call level_syncbnd_local(lev)
         return
      else if (action .eq. ELLSOLVE) then
         if (ltrace)write(*,*) 'level_apply: SendAll ELLSOLVE  ',lev
         call send_actionAll(ELLSOLVE,gi,0,lev,0,0,0,0,0,0,0)
         if (ltrace) write(*,*) 'level_apply: Calling ell_solve:'
         call ell_solve(lev)
         !do j = level_return_finest( Levelp, maxlev ), lev+1, -1
            !if(ltrace)write(*,*)'level_apply: Injecting level: ',j
            !call level_inject_local(j)
            !call level_reset_mask(j-1)
            !call level_get_resid(j-1)  ! just for debugging
         !end do
         if (ltrace) write(*,*) 'level_apply: Done w/ ELLSOLVE'
      else if (action .eq. OUTPUT) then
         if (ltrace)write(*,*) 'level_apply: SendAll OUTPUT  ',lev
         call send_actionAll(OUTPUT,gi,0,lev,0,0,0,0,0,0,0)
         ! Analysis quantities (integrals, min/max, etc):
         if (
     *      (lev.eq.0 .and. shadow.eq.0 .and. num_anafields.gt.0) .or.
     *      (lev.eq.0 .and. shadow.eq.1 .and. num_anafields.gt.0
     *                                      .and. allowedl.eq.0 ) .or.
     *      (lev.eq.1 .and. shadow.eq.1 .and. num_anafields.gt.0) ) then
            ! If on effective coarse level
            if (ltrace)write(*,*) 'level_apply: Doing analysis on: ',lev
            call level_analysis_local(lev)
         end if
         ! Surface integrals:
         if(ltrace)write(*,*)'level_apply:',asf_level,asf_period,lev
         if ( lev .eq. flev .and. asf_period.gt.0.and. .true.) then
            ! Only compute surfaces after finest level is ready for output
            tmpp   = asf_period
            ! If output requested on finer level than exists,
            ! simply output *every* step of finest level:
            if (asf_level .gt. flev) tmpp = 1
            if (asf_level .lt. flev) then
               tmpp   = asf_period*2**(flev-asf_level)
            end if
            if(ltrace)write(*,*)'level_apply:tmp:',flev,tmpp
            if(ltrace)write(*,*)'level_apply:count',lev_count(flev)
            if (mod( lev_count(flev),tmpp).eq.0) then
               if(ltrace)write(*,*)'level_apply:Computing Surface'
               !call send_actionAll(SURFCOMP,0,myid,0,0,0,0,0,0,0,0)
               call surface_compute()
               if(ltrace)write(*,*)'level_apply:Outputing Surface'
               call surface_out(grid_return_time(gi))
            end if
         end if
         !
         ! Tracer particles:  
         if (tracers_period.gt.0) then
            tmpp=tracers_period
            if ( lev .eq. 0 .and.
     .                (mod(lev_count(lev)+1,tmpp).eq.0)) then
            if(ltrace)write(*,*)'level_apply: Computing tracers'
               if(ltrace)write(*,*)'level_apply:Computing tracers'
               if(.true.)write(*,*)'level_apply:Computing tracers'
               call tracers_compute(grid_return_time(gi))
               if(ltrace)write(*,*)'level_apply:Outputing tracers'
               call tracers_out(grid_return_time(gi))
            if(ltrace)write(*,*)'level_apply: Done w/ tracers'
            end if
         end if
         !
         if(ltrace)write(*,*)'level_apply: outputting level:',lev
         call level_output_local(lev)
         return
      else if (action .eq. SETMSK) then
         ! If only one grid on level, no need to set the mask:
         if (.not.grid_return_existence(grid_return_sibling(gi))) return
         if (ltrace)write(*,*) 'level_apply: SendAll SETMSK  ',lev
         call send_actionAll(SETMSK,gi,0,lev,0,0,0,0,0,0,0)
         !
         ! Set mask of grids on level "lev" owned by this processor:
         !
         call level_set_mask(lev)
         return
      else if (action .eq. REFINE) then
         !
         ! Refinement is run by the master only:
         !
         call level_refine(lev)
         return
      else if (action .eq. INJECT) then
         if (.false.) then
            ! for debugging only:
            write(*,*) 'level_apply: Skipping INJECT on level:',lev
            return
         end if
         if (ltrace)write(*,*) 'level_apply: SendAll INJECT  ',lev
         call send_actionAll(INJECT,gi,0,lev,0,0,0,0,0,0,0)
         call level_inject_local(lev)
         call level_reset_mask(lev-1)
         return
      else if (action .eq. BOUNDS) then
         if (.false.) then
            ! for debugging only:
            write(*,*) 'level_apply: Skipping BOUNDS on level:',lev
            return
         end if
         if (ltrace)write(*,*) 'level_apply: SendAll BOUNDS  ',lev
         call send_actionAll(BOUNDS,gi,0,lev,0,0,0,0,0,0,0)
         call level_bounds_local(lev)
         if (shadow.gt.0) call level_reset_mask(lev)
         ! This sync'ing may prove necessary as a fix for
         ! a problem that develops in hyperGHEM. For now
         ! it's commented out:
         !call level_syncbnd_local(lev)
         return
      else if (action .eq. BHMASK) then
         if (ltrace)write(*,*) 'level_apply: SendAll BHMASK',lev
         call send_actionAll(BHMASK,gi,0,lev,0,0,0,0,0,0,0)
         call level_bhmask_local(lev)
         return
      else if (action .eq. INIT) then
         if (ltrace)write(*,*) 'level_apply: SendAll INIT',lev
         call send_actionAll(INIT,gi,0,lev,0,0,0,0,0,0,0)
         call level_init_local(lev)
         return
      else if (action .eq. RESETTIME) then
         if (ltrace)write(*,*) 'level_apply: SendAll RESETTIME',lev
         call send_actionAll(RESETTIME,gi,0,lev,0,0,0,0,0,0,0)
         call level_resettime_local(lev)
         return
      else if (action .eq. SAVESTATE) then
         if (ltrace)write(*,*) 'level_apply: SendAll SAVESTATE',lev
         call send_actionAll(SAVESTATE,gi,0,lev,0,0,0,0,0,0,0)
         call writestate()
         return
      else if (action .eq. READ_STATE) then
         if (ltrace)write(*,*) 'level_apply: SendAll READ_STATE',lev
         call send_actionAll(READ_STATE,gi,0,lev,0,0,0,0,0,0,0)
         call readstate()
         ! Determine some global quantities in case
         ! needed to modify data (for example, find max(rho) 
         ! to (re-)set the magnetic field:
         !call level_analysis_local(0)
         write(*,*)'level_apply: Calling level_modifydata_local()'
         call level_modifydata_local()
         write(*,*)'level_apply: Done calling level_modifydata_local()'
         return
      else if (action .eq. TRACERINIT) then
         if (ltrace)write(*,*) 'level_apply: SendAll TRACERINIT',lev
         call send_actionAll(TRACERINIT,gi,0,lev,0,0,0,0,0,0,0)
         call tracers_init()
  	  ! Variable timestep modification by Dominic Marcello 2009-08-03
      else if (action .eq. CFL) then
         if (ltrace)write(*,*) 'level_apply: SendAll CFL',lev
         call send_actionAll(CFL,gi,0,lev,0,0,0,0,0,0,0)
			call comput_cfl()
         return
	  ! --------------end modification by Dominic Marcello 2009-08-03
      else
         write(*,*) 'level_apply: *Undefined* action:',action
         write(*,*) 'level_apply: myid = ', myid
         write(*,*) 'level_apply: from = ', from
         write(*,*) 'level_apply: Quitting...'
         call my_exit('level_apply: Undefined action')
      end if

      !if (action2 .eq. SYNCBND) then
         !endtime = MPI_Wtime()
         !write(*,*) 'syncbnd took: ',endtime-starttime
      !end if


  10  format(' level_apply: lev:',I4,' grid:',I4,' action:',F4.0)
  99  format('[',I3,'] ',A,3I5)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_write_step:                                                         cc
cc                    Output a one line comment for each time step.           cc
cc                    level --- grid level that has been advanced             cc
cc                    time  --- time to which level has been advanced         cc
cc                    n     --- iteration number [0,1,2....]                  cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_write_step(level,time,n)
      implicit    none
      integer     level, n
      real(kind=8)      time
      include     'glob.inc'
      include     'grid_methods.inc'
! Variable timestep modification by Dominic Marcello 2009-08-04
      include     'param.inc'
      real        this_dt
! --------------end modification by Dominic Marcello 2009-08-04
      integer     i
      logical     ltrace
      parameter ( ltrace = .false. )
      real(kind=8)      energy,   mx, my, mz
      integer     proc_return_numprocs
      external    proc_return_numprocs

      integer        indent
      parameter    ( indent = 4 )
      character(128)  spaces
      data           spaces / ' ' /

      if (ltrace) then
         write(*,*) 'level_write_step: Enter'
         write(*,*) 'level_write_step: level = ',level
         write(*,*) 'level_write_step:  time = ',time
         write(*,*) 'level_write_step:     n = ',n
      end if


      if (level.ge.1) then

! Variable timestep modification by Dominic Marcello 2009-08-04
         if( variable_timestep .eq. 1 ) then
           this_dt = lambda * h * real(refine_factor)**(-real(level))
           write(*,1001)  spaces(1:1 + (level-1) * indent),
     &                  level, n, level_return_work(level), maxchi,
     &                  maxchi_l, time, this_dt
         else
! --------------end modification by Dominic Marcello 2009-08-04

           write(*,1000)  spaces(1:1 + (level-1) * indent),
     &                  level, n, level_return_work(level), maxchi,
     &                  maxchi_l, time

! Variable timestep modification by Dominic Marcello 2009-08-04
         endif
! --------------end modification by Dominic Marcello 2009-08-04

      else if (level.eq.0) then
         if (ltrace) write(*,*) ' Level 0 output...'
         !if (ltrace) write(*,*) ' Pointers loaded...'
c        energy = tree_integrate( 5 )  ! integrate field "edens"
c        energy = tree_integrate( 83 )  ! integrate field "edens"
         !energy = tree_integrate( 109 )  ! integrate field "edens"
         !mx     = tree_integrate( 6 )  ! integrate field "mxdens"
         !my     = tree_integrate( 7 )  ! integrate field "mydens"
         !mz     = tree_integrate( 8 )  ! integrate field "mzdens"
         energy = 0.d0
         mx     = 0.d0
         my     = 0.d0
         mz     = 0.d0
         if (ltrace) write(*,*) ' Done computing integrals'

! Variable timestep modification by Dominic Marcello 2009-08-04
         if( variable_timestep .ne. 1 ) then
! --------------end modification by Dominic Marcello 2009-08-04

           if (proc_return_numprocs().ne.1) then
           write(*,502) n,      time,     level_return_work(level),
     *                  maxchi, maxchi_l, maxchi_t
           else
           write(*,501) n,      time,     level_return_work(level),
     *                  maxchi, maxchi_l, maxchi_t,
     *                  energy, mx,       my,              mz
           if (ltrace) write(*,*) ' Done computing integDDDD'
           write(fileout,250) time,
     *                  energy, mx,       my,              mz
           end if

! Variable timestep modification by Dominic Marcello 2009-08-04
         else
           this_dt = lambda * h * real(level+1)
           if (proc_return_numprocs().ne.1) then
             write(*,602) n,      time,this_dt,level_return_work(level),
     *                    maxchi, maxchi_l, maxchi_t
           else
             write(*,601) n,      time,this_dt,level_return_work(level),
     *                    maxchi, maxchi_l, maxchi_t,
     *                    energy, mx,       my,              mz
           if (ltrace) write(*,*) ' Done computing integDDDD'
             write(fileout,250) time,
     *                    energy, mx,       my,              mz
           end if
         endif
! --------------end modification by Dominic Marcello 2009-08-04

      else 
         write(*,*) 'level_write_step: illegal level: ', level
         write(*,*) 'level_write_step: Quitting...'
         Stop
      end if

 1000 format('     ',A,'|',I2,'| (',I6,') W:',F8.0,' [',G8.1,'/',I2,
     *             '] t:',F9.3)
! Variable timestep modification by Dominic Marcello 2009-08-04
 1001 format('     ',A,'|',I2,'| (',I4,') W:',F8.0,' [',G8.1,'/',I2,
     *             '] t:',ES10.3, ' dt:',ES10.3)
 601  format(' Step ',i5,' t: ',ES10.3,' dt: ', es10.3, 
     *       ' W:',F5.0,' MaxChi:',G8.1,'/',I2,'/',F5.2,' E: ',4f6.1)
 602  format(' Step ',i5,' t: ',es10.3, ' dt: ', es10.3, 
     *       ' W:',F8.0,' MaxChi:',G8.1,'/',I2,'/',F5.2)
 ! --------------end modification by Dominic Marcello 2009-08-04
 501  format(' Step ',i5,' t: ',f8.3,' W:',F8.0,' MaxChi:',G8.1,'/',
     *       I2,'/',F5.2,' E: ',4f6.1)
 502  format(' Step ',i5,' t: ',f8.3,' W:',F8.0,' MaxChi:',G8.1,'/',
     *       I2,'/',F5.2)
 250  format(E15.7,
     *       ' E: ',E24.16,' MX: ',E24.16,' MY: ',E24.16,' MZ: ',E24.16)

      if (ltrace) then
         write(*,*) 'level_write_step: Exit'
      end if

      return
      end     ! END: level_write_step

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_find_bounds:                                                        cc
cc                    Get min/max bounds for an entire level.                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_find_bounds(level,minx,maxx,miny,maxy,minz,maxz)
      implicit     none
      integer      level
      real(kind=8) minx,maxx, miny,maxy, minz,maxz
      include     'grid_methods.inc'
      real(kind=8) tminx,tmaxx, tminy,tmaxy, tminz,tmaxz
      integer      gi

      minx = 0.d0
      miny = 0.d0
      minz = 0.d0
      maxx = 0.d0
      maxy = 0.d0
      maxz = 0.d0

      gi = level_return_start(level)
      if (.not.grid_return_existence(gi)) then
         write(*,*)'level_find_bounds:Problem: Level doesnt exist',level
         return
      end if
      call grid_find_bounds(gi,minx,maxx,miny,maxy,minz,maxz)

 10   gi = grid_return_sibling(gi)
      if (grid_return_existence(gi)) then
         call grid_find_bounds(gi,tminx,tmaxx,tminy,tmaxy,tminz,tmaxz)
         if (tminx .lt. minx) minx = tminx
         if (tminy .lt. miny) miny = tminy
         if (tminz .lt. minz) minz = tminz
         if (tmaxx .gt. maxx) maxx = tmaxx
         if (tmaxy .gt. maxy) maxy = tmaxy
         if (tmaxz .gt. maxz) maxz = tmaxz
         goto 10
      end if

      return 
      end    ! END: level_find_bounds

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_mg_cgccompute_local:                                                cc
cc                    Call grid_relax() on all local grids of a level.        cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_mg_cgccompute_local(level)
      implicit     none
      integer      level
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      gi

      gi = level_return_start(level)
      !
  10  if (grid_return_existence(gi)) then
         if (grid_is_local(gi)) then
            call grid_cgccompute(gi)
         end if
         gi = grid_return_sibling(gi)
         goto 10
      end if

      return 
      end    ! END: level_mg_cgccompute_local

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_mg_pstcompute_local:                                                cc
cc                    Call grid_relax() on all local grids of a level.        cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_mg_pstcompute_local(level)
      implicit     none
      integer      level
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      gi

      gi = level_return_start(level)
      !
  10  if (grid_return_existence(gi)) then
         if (grid_is_local(gi)) then
            call grid_pstcompute(gi)
         end if
         gi = grid_return_sibling(gi)
         goto 10
      end if

      return 
      end    ! END: level_mg_pstcompute_local

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_mgsmooth:                                                           cc
cc                  Call grid_mgsmooth()on all local grids of a level.        cc
cc                  This is an attempt to make things more robust where       cc
cc                  when nonverging of vcycles is detected, the guess is      cc
cc                  smoothed.                                                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_mgsmooth(level,numattempts)
      implicit     none
      integer      level, numattempts
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      gi

      gi = level_return_start(level)
      !
  10  if (grid_return_existence(gi)) then
         if (grid_is_local(gi)) then
            call grid_mgsmooth(gi, numattempts)
         end if
         gi = grid_return_sibling(gi)
         goto 10
      end if

      return 
      end    ! END: level_mgsmooth

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_relax_local:                                                        cc
cc                    Call grid_relax() on all local grids of a level.        cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_relax_local(level)
      implicit     none
      integer      level
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      gi

      gi = level_return_start(level)
      !
  10  if (grid_return_existence(gi)) then
         if (grid_is_local(gi)) then
            call grid_relax(gi)
         end if
         gi = grid_return_sibling(gi)
         goto 10
      end if

      return 
      end    ! END: level_relax_local

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_iter_local:                                                         cc
cc                    Call grid_iter() all local grids of a level.            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_iter_local(level)
      implicit     none
      integer      level
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      include     'glob.inc'
      include     'grid.inc'
      integer      gi
      ! Field number to be output if tracing is on:
      character(2) iterch
      integer     outfnum
      parameter ( outfnum  =  1 )
      !parameter ( outfnum  = 51 )
      real(kind=8) orignorm,newnorm
      real(kind=8) myl2norm3d
      external     myl2norm3d
      logical      ltrace
      parameter (  ltrace = .false.) 
      logical      ltrace2
      parameter (  ltrace2 = .false.) 

      !ltrace = .false.
      !if (level.eq.6) ltrace = .true.

      if (ltrace) then
         write(*,99)myid,'level_iter: Enter on level:',level
         write(*,99)myid,'level_iter: Run grid_test on grids',level
         call level_test_local()
      end if

      gi = level_return_start(level)
      !
  10  if (grid_return_existence(gi)) then
         if(ltrace)write(*,99)myid,'level_iter_local:gi',gi
         if (grid_is_local(gi)) then
          if(ltrace)write(*,99)myid,'level_iter_local:gi',gi,gr_iter(gi)
            orignorm = myl2norm3d(q(gfunc_pointer(outfnum+1,gi)),
     *                                 gr_nx(gi),gr_ny(gi),gr_nz(gi))
            if (ltrace2) then
            !if (gi.eq.19.and.ltrace) then
               call int2str(gr_iter(gi),iterch)
               call field_out3d(q(gfunc_pointer(outfnum,gi)),
     *             gr_t(gi),'outfnum'//iterch,
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             gr_nx(gi),gr_ny(gi),gr_nz(gi),myid)
               call field_out3d(q(gfunc_pointer(outfnum+1,gi)),
     *             gr_t(gi),'outfnumADV'//iterch,
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             gr_nx(gi),gr_ny(gi),gr_nz(gi),myid)
            end if
            if (ltrace) then
              write(*,99)myid,'level_iter:XRun grid_test on:',level,gi
              call level_test_local()
            end if
            if(ltrace2)call grid_dump_info(gi)
            if(ltrace)write(*,99)myid,'level_iter:calling grid_iter'
            call grid_iter(gi)
            if(ltrace)write(*,99)myid,'level_iter:Back from grid_iter'
            if(ltrace2)call grid_dump_info(gi)
            if (ltrace) then
              write(*,99)myid,'level_iter:YRun grid_test on:',level,gi
              call level_test_local()
            end if
            if (assume_symmetry.ne.0) then
               if(ltrace)write(*,99)myid,'level_iter_local:grid_symm',gi
               call grid_symm(gi)
            end if
            if (ltrace2) then
               newnorm=myl2norm3d(q(gfunc_pointer(outfnum+1,gi)),
     *                                 gr_nx(gi),gr_ny(gi),gr_nz(gi))
               if (abs(newnorm-orignorm).gt.newnorm*100.0) then
                 write(*,99)myid,'level_iter_local: PROBLEM w/gi: ',gi
               write(*,99)myid,'level_iter_local: outfnum+1: ',outfnum+1
                 write(*,98)myid,'level_iter_local: orig/newnorm: ',
     *                                 orignorm,newnorm
                 write(*,99)myid,'level_iter_local:gi',gi,gr_iter(gi)
                 write(*,99)myid,'level_iter_local:Call grid_dump_info'
                 call grid_dump_info(gi)
               end if
               call field_out3d(q(gfunc_pointer(outfnum,gi)),
     *             gr_t(gi),'outfnumPOS'//iterch,
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             gr_nx(gi),gr_ny(gi),gr_nz(gi),myid)
               call field_out3d(q(gfunc_pointer(outfnum+1,gi)),
     *             gr_t(gi),'outfnumADVPOS'//iterch,
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             gr_nx(gi),gr_ny(gi),gr_nz(gi),myid)
            end if
         else
            if(ltrace) write(*,99)myid,'level_iter:Gridnotlocal',gi
         end if   ! if grid is local
         if (ltrace) then
           write(*,99)myid,'level_iter:ZRun grid_test on grids',level,gi
           call level_test_local()
         end if
         gi = grid_return_sibling(gi)
         goto 10
      end if

      if (ltrace) then
         write(*,99)myid,'level_iter:last grid_test on grids',level
         call level_test_local()
         write(*,99)myid,'level_iter: Done on level:',level
      end if

  98  format('[',I4,'] ',A,3G15.7)
  99  format('[',I4,'] ',A,3I5)
      return 
      end    ! END: level_iter_local

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_step2_local:                                                        cc
cc                    Call grid_step2() all local grids of a level.           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_step2_local(level)
      implicit     none
      integer      level
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      gi

      gi = level_return_start(level)
      !
  10  if (grid_return_existence(gi)) then
         if (grid_is_local(gi)) then
            call grid_step2(gi)
         end if
         gi = grid_return_sibling(gi)
         goto 10
      end if

      return 
      end    ! END: level_step2_local

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_reset_mask:                                                         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_reset_mask(level)
      implicit     none
      integer      level
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      gi

      gi = level_return_start(level)
      !
  10  if (grid_return_existence(gi)) then
         if (grid_is_local(gi)) then
            call grid_reset_mask(gi)
         end if
         gi = grid_return_sibling(gi)
         goto 10
      end if

      return 
      end    ! END: level_reset_mask

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_ellpost_local:                                                      cc
cc                       Have each grid on levels as fine or finer than       cc
cc                       than "lev" postprocess what it needs to              cc
cc                                   after the elliptic solver completes.     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_ellpost_local(lev)
      implicit     none
      integer      lev
      include     'glob.inc'
      include     'grid_methods.inc'
      integer      gi
      integer      li
      logical      ltrace
      parameter (  ltrace = .false. )

      if (ltrace) write(*,*) 'level_ellpost_local: Finding post:'
      li = lev
      !li = 0
      do while (Levelp(li).gt.0)
         if (ltrace) write(*,*) 'level_ellpost_local: li= ',li
         gi = level_return_start(li)
         !
  10     if (grid_return_existence(gi)) then
            if (ltrace) write(*,*) 'level_ellpost_local: gi= ',gi
            if (grid_is_local(gi)) then
               if (ltrace) write(*,*) 'level_ellpost_local: Local'
               call grid_ellpost(gi)
               call grid_comp_deriv(gi)
            end if
            gi = grid_return_sibling(gi)
            goto 10
         end if
         li = li + 1
      end do

      return 
      end    ! END: level_cellpost_local

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_ellrhs_local:                                                       cc
cc                       Have each grid on a   levels compute the rhs for     cc
cc                       the elliptic solver.                                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_ellrhs_local(lev)
      implicit     none
      integer      lev
      include     'glob.inc'
      include     'grid_methods.inc'
      integer      gi
      !integer      li
      logical      ltrace
      parameter (  ltrace = .false. )

      if (ltrace) write(*,*) 'level_ellrhs_local: Finding rhs:'
      !li = lev
      !do while (Levelp(li).gt.0)
         !if (ltrace) write(*,*) 'level_ellrhs_local: li= ',li
         !gi = level_return_start(li)
         gi = level_return_start(lev)
         !
  10     if (grid_return_existence(gi)) then
            if (ltrace) write(*,*) 'level_ellrhs_local: gi= ',gi
            if (grid_is_local(gi)) then
               if (ltrace) write(*,*) 'level_ellrhs_local: Local'
               call grid_ellrhs(gi)
            end if
            gi = grid_return_sibling(gi)
            goto 10
         end if
         !li = li + 1
      !end do

      return 
      end    ! END: level_ellrhs_local

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_comput_local:                                                       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_comput_local(level)
      implicit     none
      integer      level
      !include     'mpif.h'
      !include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      gi

      gi = level_return_start(level)
      !
  10  if (grid_return_existence(gi)) then
         if (grid_is_local(gi)) then
            call grid_comp_deriv(gi)
         end if
         gi = grid_return_sibling(gi)
         goto 10
      end if

      return 
      end    ! END: level_comput_local

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_output_local:                                                       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_output_local(level)
      implicit     none
      integer      level
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      gi

      logical      ltrace
      parameter (  ltrace = .false. )

      if (ltrace) write(*,99) myid, 'level_output_local: level= ',level

      gi = level_return_start(level)
      !
  10  if (grid_return_existence(gi)) then
         if (grid_is_local(gi)) then
            if (ltrace) write(*,99) myid, 'level_output_local: gi=',gi
            call grid_output_sdf(gi)
         end if
         gi = grid_return_sibling(gi)
         goto 10
      end if

      if (ltrace) write(*,99) myid, 'level_output_local: Finished'

  99  format('[',I3,'] ',A,3I5)
      return 
      end    ! END: level_output_local

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_syncbnd_local:                                                      cc
cc                    Loop over level and syncbnd all grids local to this proccc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_syncbnd_local(level)
      implicit     none
      integer      level
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      gi, gj, gi_owner, gj_owner
      integer      round, numrounds, proc
      integer      round_robin
      external     round_robin
      logical      gi_local, gj_local
      logical      ltrace
      parameter (  ltrace = .false. )
      ! Use large messages in sync?
      logical      synclarge
      parameter (  synclarge = .false. )

      if (.false.) then
         write(*,99) myid, 'level_syncbnd_local: Turning off.'
         return
      end if

      if (ltrace) then
         write(*,99) myid, 'level_syncbnd_local: level = ',level
         write(*,99)myid,'level_syncbnd: Run level_test_local'
         call level_test_local()
      end if

      ! Used to test the whole round_robin routine:
      !call test_round_robin(7)

      !
      ! We are processor "myid" and we need to check with all other
      ! processors to see if we need to sync with them. We want to do
      ! so in a way that avoids having a processor wait for the other
      ! processor to finish with someone else. So we use a round robin.
      ! 
      ! For each "round" of the round robin, we determine which processor
      ! with which to sync.

      if (mod(numprocs,2).eq.0) then
         ! EVEN number of processors:
         numrounds = numprocs - 1
      else
         ! ODD  number of processors:
         numrounds = numprocs
      end if
      ! Add a round to sync pairs of grids both owned by this process:
      numrounds = numrounds + 1

      !
      ! Loop over the rounds of the round robin:
      !
      do round = 1, numrounds
        !
        ! Determine which processor with which to sync:
        !
        if (round .ne. numrounds) then
           proc = round_robin(round, myid, numprocs)
        else
           proc = myid
        end if
        !
        ! Make sure we don't have a "bye" (only happens for odd numprocs):
        !
        if (proc .ge. 0) then
           if (ltrace)
     *        write(*,99) myid, 'level_syncbnd_local: Sync w/proc=',proc
           !
           ! Loop for all grids (gi) on level. For each gi,
           ! consider pairings of gi and gj where occurs in the
           ! hierarchy after gi (this avoids duplicate pairings),
           ! and sync those pairings (gi,gj):
           !
           gi = level_return_start(level)
           !
 10        if (grid_return_existence(gi)) then
              if (ltrace) write(*,99) myid, '      gi = ',gi
              gi_owner = grid_return_owner(gi)
              gi_local = grid_is_local(gi)
              ! If neither us nor processor "proc" own gi, then skip:
              if (.not. gi_local .and. gi_owner .ne. proc) goto 40
              !
              gj = grid_return_sibling(gi)
 20           if (grid_return_existence(gj)) then
                 !
                 ! Should we skip this pairing?
                 gj_owner = grid_return_owner(gj)
                 gj_local = grid_is_local(gj)
                 if (.not. gj_local .and. gj_owner .ne. proc) goto 30
                 if (gi_local.and.gj_local .and. (myid.ne.proc)) goto 30
                 ! Does pair intersect?
                 if (ltrace) write(*,99) myid, '      w/ gj = ',gj
                 if ( grid_intersect(gi,gj) ) then
                    if (gj_local) then
                       if (gi_local) then
                          if (ltrace) write(*,99) myid, '   gi/gj local'
                          call grid_syncbnd(gj, gi)
                       else
                          if (ltrace) write(*,99) myid, '   gj   local'
                          if (synclarge) then
                             call grid_syncbnd_slaveLG(gj, gi)
                          else
                             call grid_syncbnd_slave(gj, gi)
                          endif
                       end if
                    else
                       if (gi_local) then
                          if (ltrace) write(*,99) myid, '   gi   local'
                          if (synclarge) then
                             call grid_syncbnd_slaveLG(gi, gj)
                          else
                             call grid_syncbnd_slave(gi, gj)
                          endif
                       else
                          if (ltrace) write(*,99) myid, '   none local'
                       end if
                    end if
                 else
                          if (ltrace) write(*,99) myid, '   No intersec'
                 end if
                 ! Look to next sibling
  30             gj = grid_return_sibling(gj)
                 goto 20
              end if
              !
              ! No more siblings w/ which to sync gi,
              ! so, repeat process for next grid on level:
              !
  40          gi = grid_return_sibling(gi)
              goto 10
           end if

        end if
      end do

      if (ltrace) then
         write(*,99)myid,'level_syncbnd: Run level_test_local'
         call level_test_local()
      end if

      if (ltrace) write(*,99) myid, 'level_syncbnd_loca: Done.'
  99  format('[',I4,'] ',A,3I5)

      return 
      end    ! END: level_syncbnd_local

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_syncbndi_local:                                                     cc
cc                    Sync bounds *just* for *elliptic* fields.               cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_syncbndi_local(level)
      implicit     none
      integer      level
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      gi, gj, gi_owner, gj_owner
      integer      round, numrounds, proc
      integer      round_robin
      external     round_robin
      logical      gi_local, gj_local
      logical      ltrace
      parameter (  ltrace = .false. )

      if (.false.) then
         write(*,99) myid, 'level_syncbndi_local: Turning off.'
         return
      end if

      if (ltrace) then
         write(*,99) myid, 'level_syncbndi_local: level = ',level
      end if

      ! Used to test the whole round_robin routine:
      !call test_round_robin(7)

      !
      ! We are processor "myid" and we need to check with all other
      ! processors to see if we need to sync with them. We want to do
      ! so in a way that avoids having a processor wait for the other
      ! processor to finish with someone else. So we use a round robin.
      ! 
      ! For each "round" of the round robin, we determine which processor
      ! with which to sync.

      if (mod(numprocs,2).eq.0) then
         ! EVEN number of processors:
         numrounds = numprocs - 1
      else
         ! ODD  number of processors:
         numrounds = numprocs
      end if
      ! Add a round to sync pairs of grids both owned by this process:
      numrounds = numrounds + 1

      !
      ! Loop over the rounds of the round robin:
      !
      do round = 1, numrounds
        !
        ! Determine which processor with which to sync:
        !
        if (round .ne. numrounds) then
           proc = round_robin(round, myid, numprocs)
        else
           proc = myid
        end if
        !
        ! Make sure we don't have a "bye" (only happens for odd numprocs):
        !
        if (proc .ge. 0) then
           if (ltrace)
     *        write(*,99) myid,'level_syncbndi_local: Sync w/proc=',proc
           !
           ! Loop for all grids (gi) on level. For each gi,
           ! consider pairings of gi and gj where occurs in the
           ! hierarchy after gi (this avoids duplicate pairings),
           ! and sync those pairings (gi,gj):
           !
           gi = level_return_start(level)
           !
 10        if (grid_return_existence(gi)) then
              if (ltrace) write(*,99) myid, '      gi = ',gi
              gi_owner = grid_return_owner(gi)
              gi_local = grid_is_local(gi)
              ! If neither us nor processor "proc" own gi, then skip:
              if (.not. gi_local .and. gi_owner .ne. proc) goto 40
              !
              gj = grid_return_sibling(gi)
 20           if (grid_return_existence(gj)) then
                 !
                 ! Should we skip this pairing?
                 gj_owner = grid_return_owner(gj)
                 gj_local = grid_is_local(gj)
                 if (.not. gj_local .and. gj_owner .ne. proc) goto 30
                 if (gi_local.and.gj_local .and. (myid.ne.proc)) goto 30
                 ! Does pair intersect?
                 if (ltrace) write(*,99) myid, '      w/ gj = ',gj
                 if ( grid_intersect(gi,gj) ) then
                    if (gj_local) then
                       if (gi_local) then
                          if (ltrace) write(*,99) myid, '   gi/gj local'
                          call grid_syncbndi(gj, gi)
                       else
                          if (ltrace) write(*,99) myid, '   gj   local'
                          call grid_syncbndi_slave(gj, gi)
                       end if
                    else
                       if (gi_local) then
                          if (ltrace) write(*,99) myid, '   gi   local'
                          call grid_syncbndi_slave(gi, gj)
                       else
                          if (ltrace) write(*,99) myid, '   none local'
                       end if
                    end if
                 else
                          if (ltrace) write(*,99) myid, '   No intersec'
                 end if
                 ! Look to next sibling
  30             gj = grid_return_sibling(gj)
                 goto 20
              end if
              !
              ! No more siblings w/ which to sync gi,
              ! so, repeat process for next grid on level:
              !
  40          gi = grid_return_sibling(gi)
              goto 10
           end if

        end if
      end do

      if (ltrace) write(*,99) myid, 'level_syncbndi_loca: Done.'
  99  format('[',I3,'] ',A,3I5)

      return 
      end    ! END: level_syncbndi_local

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_mginject_local:                                                     cc
cc                    Loop over level and inject all grids local to this proc cc
cc               NB:  copied from level_inject_local()                        cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_mginject_local(level)
      implicit     none
      integer      level
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      gi, gj, gi_owner, gj_owner
      integer      round, numrounds, proc
      integer      round_robin
      external     round_robin
      logical      gi_local, gj_local
      logical      ltrace
      parameter (  ltrace = .false. )

      if (.false.) then
         write(*,99) myid, 'level_mginject_local: Turning off.'
         return
      end if

      if (ltrace) then
         write(*,99) myid, 'level_mginject_local: level = ',level
      end if

      !
      ! We are processor "myid" and we need to check with all other
      ! processors to see if we need to sync with them. We want to do
      ! so in a way that avoids having a processor wait for the other
      ! processor to finish with someone else. So we use a round robin.
      ! 
      ! For each "round" of the round robin, we determine which processor
      ! with which to sync.

      if (mod(numprocs,2).eq.0) then
         ! EVEN number of processors:
         numrounds = numprocs - 1
      else
         ! ODD  number of processors:
         numrounds = numprocs
      end if
      ! Add a round to sync pairs of grids both owned by this process:
      numrounds = numrounds + 1

      !
      ! Loop over the rounds of the round robin:
      !
      do round = 1, numrounds
        if(ltrace)write(*,99) myid, 'level_mginject_local:round=',round
        !
        ! Determine which processor with which to sync:
        !
        if (round .ne. numrounds) then
           proc = round_robin(round, myid, numprocs)
        else
           proc = myid
        end if
        !
        ! Make sure we don't have a "bye" (only happens for odd numprocs):
        !
        if (proc .ge. 0) then
           if (ltrace) then
              write(*,99)myid,'level_mginject_local:Inject w/proc=',proc
           end if
           !
           ! Loop for all grids (gi) on level. For each gi,
           ! consider pairings of gi and gj where occurs in the
           ! hierarchy after gi (this avoids duplicate pairings),
           ! and sync those pairings (gi,gj):
           !
           gi = level_return_start(level)
           !
 10        if (grid_return_existence(gi)) then
              if (ltrace) write(*,99) myid, '      gi = ',gi
              gi_owner = grid_return_owner(gi)
              gi_local = grid_is_local(gi)
              ! If neither us nor processor "proc" own gi, then skip:
              if (.not. gi_local .and. gi_owner .ne. proc) goto 40
              !
              !gj = grid_return_sibling(gi)
              ! Loop over parent level:
              gj = level_return_start(level-1)
 20           if (grid_return_existence(gj)) then
                 !
                 ! Should we skip this pairing?
                 gj_owner = grid_return_owner(gj)
                 gj_local = grid_is_local(gj)
                 if (.not. gj_local .and. gj_owner .ne. proc) goto 30
                 if (gi_local.and.gj_local .and. (myid.ne.proc)) goto 30
                 ! Does pair intersect?
                 if (ltrace) write(*,99) myid, '      w/ gj = ',gj
                 if ( grid_intersect(gi,gj) ) then
                    if (gj_local) then
                       if (gi_local) then
                          if (ltrace) write(*,99) myid, '   gi/gj local'
                          call grid_mginject(gi, gj)
                       else
                          if (ltrace) write(*,99) myid, '   gj   local'
                          call grid_mginject_master(gi, gj)
                       end if
                    else
                       if (gi_local) then
                          if (ltrace) write(*,99) myid, '   gi   local'
                          call grid_mginject_slave(gi, gj)
                       else
                          if (ltrace) write(*,99) myid, '   none local'
                       end if
                    end if
                 else
                          if (ltrace) write(*,99) myid, '   No intersec'
                 end if
                 ! Look to next sibling
  30             gj = grid_return_sibling(gj)
                 goto 20
              end if
              !
              ! No more siblings w/ which to sync gi,
              ! so, repeat process for next grid on level:
              !
  40          gi = grid_return_sibling(gi)
              goto 10
           end if

        end if
      end do

  99  format('[',I3,'] ',A,3I5)

      return 
      end    ! END: level_mginject_local

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_inject_local:                                                       cc
cc                    Loop over level and inject all grids local to this proc cc
cc               NB:  similar to level_syncbnd_local()                        cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_inject_local(level)
      implicit     none
      integer      level
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      gi, gj, gi_owner, gj_owner
      integer      round, numrounds, proc
      integer      round_robin
      external     round_robin
      logical      gi_local, gj_local
      logical      ltrace
      parameter (  ltrace = .false. )

      if (.false.) then
         write(*,99) myid, 'level_inject_local: Turning off.'
         return
      end if

      if (ltrace) then
         write(*,99) myid, 'level_inject_local: level = ',level
      end if

      !
      ! We are processor "myid" and we need to check with all other
      ! processors to see if we need to sync with them. We want to do
      ! so in a way that avoids having a processor wait for the other
      ! processor to finish with someone else. So we use a round robin.
      ! 
      ! For each "round" of the round robin, we determine which processor
      ! with which to sync.

      if (mod(numprocs,2).eq.0) then
         ! EVEN number of processors:
         numrounds = numprocs - 1
      else
         ! ODD  number of processors:
         numrounds = numprocs
      end if
      ! Add a round to sync pairs of grids both owned by this process:
      numrounds = numrounds + 1

      !
      ! Loop over the rounds of the round robin:
      !
      do round = 1, numrounds
        if (ltrace) write(*,99) myid, 'level_inject_local:round=',round
        !
        ! Determine which processor with which to sync:
        !
        if (round .ne. numrounds) then
           proc = round_robin(round, myid, numprocs)
        else
           proc = myid
        end if
        !
        ! Make sure we don't have a "bye" (only happens for odd numprocs):
        !
        if (proc .ge. 0) then
           if (ltrace) then
              write(*,99) myid, 'level_inject_local:Inject w/proc=',proc
           end if
           !
           ! Loop for all grids (gi) on level. For each gi,
           ! consider pairings of gi and gj where occurs in the
           ! hierarchy after gi (this avoids duplicate pairings),
           ! and sync those pairings (gi,gj):
           !
           gi = level_return_start(level)
           !
 10        if (grid_return_existence(gi)) then
              if (ltrace) write(*,99) myid, '      gi = ',gi
              gi_owner = grid_return_owner(gi)
              gi_local = grid_is_local(gi)
              ! If neither us nor processor "proc" own gi, then skip:
              if (.not. gi_local .and. gi_owner .ne. proc) goto 40
              !
              !gj = grid_return_sibling(gi)
              ! Loop over parent level:
              gj = level_return_start(level-1)
 20           if (grid_return_existence(gj)) then
                 !
                 ! Should we skip this pairing?
                 gj_owner = grid_return_owner(gj)
                 gj_local = grid_is_local(gj)
                 if (.not. gj_local .and. gj_owner .ne. proc) goto 30
                 if (gi_local.and.gj_local .and. (myid.ne.proc)) goto 30
                 ! Does pair intersect?
                 if (ltrace) write(*,99) myid, '      w/ gj = ',gj
                 if ( grid_intersect(gi,gj) ) then
                    if (gj_local) then
                       if (gi_local) then
                          if (ltrace) write(*,99) myid, '   gi/gj local'
                          call grid_inject(gi, gj)
                       else
                          if (ltrace) write(*,99) myid, '   gj   local'
                          call grid_inject_master(gi, gj)
                       end if
                    else
                       if (gi_local) then
                          if (ltrace) write(*,99) myid, '   gi   local'
                          call grid_inject_slave(gi, gj)
                       else
                          if (ltrace) write(*,99) myid, '   none local'
                       end if
                    end if
                 else
                          if (ltrace) write(*,99) myid, '   No intersec'
                 end if
                 ! Look to next sibling
  30             gj = grid_return_sibling(gj)
                 goto 20
              end if
              !
              ! No more siblings w/ which to sync gi,
              ! so, repeat process for next grid on level:
              !
  40          gi = grid_return_sibling(gi)
              goto 10
           end if

        end if
      end do

      if (ltrace) then
         write(*,99)myid,'level_syncbnd: Run grid_test on grids',level
         call level_test_local()
      end if

  99  format('[',I3,'] ',A,3I5)

      return 
      end    ! END: level_inject_local

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_bounds_local:                                                       cc
cc                    Loop over level and get bounds all grids local to proc  cc
cc               NB:  similar to level_inject_local()                         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_bounds_local(level)
      implicit     none
      integer      level
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      include     'glob.inc'          ! for shadow
      integer      gi, gj, gi_owner, gj_owner
      integer      round, numrounds, proc
      integer      round_robin
      external     round_robin
      logical      gi_local, gj_local
      logical      first_time_through
      logical      ltrace
      parameter (  ltrace  = .false.  )
      logical      ltrace2
      parameter (  ltrace2 = .false. )
      logical      debug
      parameter (  debug   = .false. )

      if (.false.) then
         write(*,99) myid, 'level_bounds_local: Turning off.'
         return
      end if

      if (ltrace2) then
         write(*,99) myid, 'level_bounds_local: level = ',level
      end if

      if (debug) then
         write(*,99)myid,'level_bounds_local: Setting points'
         write(*,99)myid,'level_bounds_local:to be re-set toLARGENUMBER'
         write(*,99)myid,'level_bounds_local:afterwards,should be reset'
         gi = level_return_start(level)
         !
 90      if (grid_return_existence(gi)) then
            if(ltrace)write(*,99)myid,'level_bounds_local:gi:',gi
            if (grid_is_local(gi)) then
               if(.true.)write(*,99)myid,'level_bounds_local:local' 
               call grid_setbndslarge(gi)
            end if
            gi = grid_return_sibling(gi)
            goto 90
         end if
         if(ltrace)write(*,99)myid,'level_bounds_local:Done w/LARGE'
      end if

      !
      ! We are processor "myid" and we need to check with all other
      ! processors to see if we need to sync with them. We want to do
      ! so in a way that avoids having a processor wait for the other
      ! processor to finish with someone else. So we use a round robin.
      ! 
      ! For each "round" of the round robin, we determine which processor
      ! with which to sync.

      if (mod(numprocs,2).eq.0) then
         ! EVEN number of processors:
         if (ltrace) write(*,95)myid, '   Numprocs EVEN: ',numprocs
         numrounds = numprocs - 1
      else
         ! ODD  number of processors:
         if (ltrace) write(*,95)myid, '   Numprocs ODD:  ',numprocs
         numrounds = numprocs
      end if
      ! Add a round to sync pairs of grids both owned by this process:
      numrounds = numrounds + 1

      if (ltrace2) then
         write(*,99) myid, 'level_bounds_local: level     = ',level
         write(*,99) myid, 'level_bounds_local: numprocs  = ',numprocs
         write(*,99) myid, 'level_bounds_local: numrounds = ',numrounds
      end if

      !
      ! Loop over the rounds of the round robin:
      !
      first_time_through = .true.
      do round = 1, numrounds
        !
        ! Determine which processor with which to sync:
        !
        if (round .ne. numrounds) then
           proc = round_robin(round, myid, numprocs)
        else
           proc = myid
        end if
        !
        ! Make sure we don't have a "bye" (only happens for odd numprocs):
        !
        if (proc .ge. 0) then
           if (ltrace2)
     *        write(*,99) myid, 'level_bounds_local:Bounds w/proc=',proc
           !
           ! Loop for all grids (gi) on level. For each gi,
           ! consider pairings of gi and gj where occurs in the
           ! hierarchy after gi (this avoids duplicate pairings),
           ! and sync those pairings (gi,gj):
           !
           gi = level_return_start(level)
           !
 10        if (grid_return_existence(gi)) then
              if (ltrace2) write(*,99) myid, '      gi = ',gi
              gi_owner = grid_return_owner(gi)
              gi_local = grid_is_local(gi)
              ! If neither us nor processor "proc" own gi, then skip:
              if (.not. gi_local .and. gi_owner .ne. proc) goto 40
              if ( first_time_through .and.
     *             gi_local           .and.
     *             shadow.gt.0              ) then
                 if(ltrace)write(*,99)myid,'    Zero error on grid:',gi
                 call grid_zeroerror(gi)
              end if
              !
              ! Loop over parent level:
              gj = level_return_start(level-1)
 20           if (grid_return_existence(gj)) then
                 !
                 ! Should we skip this pairing?
                 gj_owner = grid_return_owner(gj)
                 gj_local = grid_is_local(gj)
                 ! Look at next sibling if neither proc owns it:
                 if (.not. gj_local .and. gj_owner .ne. proc) goto 30
                 ! Both grids owned by this processor, but don't 
                 ! carryout unless we're in that round:
                 if (gi_local.and.gj_local .and. (myid.ne.proc)) goto 30
                 ! Does pair intersect?
                 if (ltrace2)write(*,99) myid, '      gi w/ gj = ',gi,gj
                 if (ltrace) write(*,95)myid, 'round: ',round,' proc: ',
     *               proc,' gi: ',gi, ' gj: ',gj
                 if ( grid_intersect(gi,gj) ) then
                    if (gj_local) then
                       if (gi_local) then
                          if (ltrace) write(*,99) myid, '   gi/gj local'
                          call grid_intp_bnd(gi, gj)
                       else
                          if (ltrace) write(*,99) myid, '   gj   local'
                          call grid_intp_bnd_master(gi, gj)
                       end if
                    else
                       if (gi_local) then
                          if (ltrace) write(*,99) myid, '   gi   local'
                          call grid_intp_bnd_slave(gi, gj)
                       else
                          if (ltrace) write(*,99) myid, '   none local'
                       end if
                    end if
                 else
                          if (ltrace) write(*,99) myid, '  No intersec'
                 end if
                 ! Look to next sibling
  30             gj = grid_return_sibling(gj)
                 goto 20
              end if
              !
              ! No more siblings w/ which to sync gi,
              ! so, repeat process for next grid on level:
              !
  40          gi = grid_return_sibling(gi)
              goto 10
           end if
           !
           first_time_through = .false.
        end if
      end do

      if (ltrace) then
         write(*,99)myid,'level_bounds: Run grid_test on grids',level
         call level_test_local()
      end if

      if (ltrace2) then
         write(*,99) myid, 'level_bounds_local: Finished'
      end if

  95  format('[',I3,'] ',A,I5,A,I5,A,I5,A,I5)
  99  format('[',I3,'] ',A,3I5)

      return 
      end    ! END: level_bounds_local

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_bhmask_local:                                                       cc
cc                       Each processor calls this at the initial time        cc
cc                  and after each time step of the finest level in           cc
cc                  existence with the assumption that the masked regions     cc
cc                  are completely contained within this finest level.        cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_bhmask_local(level)
      implicit     none
      integer      level
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      include     'glob.inc'
      include     'grid.inc'    ! needed for debugging output only
      include     'mask.inc'
      include     'param.inc'
      integer      gi, i, k, kk, ii
      integer      tmp_num_masks
      logical      tmp_bh_true(max_num_masks)
      real(kind=8) time, hi, motion,motionx,motiony, omega
      real(kind=8) distance, mindistance
      real(kind=8) rad_left_x,rad_left_y,rad_left_z       
      real(kind=8) rad_right_x,rad_right_y,rad_right_z        
      real(kind=8) tmppoint(3,max_num_masks), loc_min(max_num_masks),
     *                                       glob_min(max_num_masks)
      integer      proc_min(max_num_masks)
      character(128) fname,  procname
      integer      mylen
      ! Simple level counter when looping over levels:
      integer      lev_i,finestlev
      real(kind=8) oldmincoord, oldmaxcoord
      real(kind=8) newmincoord, newmaxcoord, newradius
      logical      noprevioushole(max_num_masks)

      logical      mask_changed, any_mask_changed, levelchanged
      logical      double_equal
      external     double_equal
      include     'largesmall.inc'
      real(kind=8) small_delta, span, weight

      !
      ! By how much to shrink "horizon" volume to get excision volume:
      !                ( 0 < SHRINK < 1 )
      !
      real(kind=8) SHRINK
      !parameter (  SHRINK = 0.30d0 )
      !parameter (  SHRINK = 0.20d0 )
      parameter (  SHRINK = 0.8d0 )
      !
      ! Maximum number of times to iterate to find
      ! the right center for the hole:
      !
      integer      iter
      integer      MAXITERS
      parameter (  MAXITERS = 1 )
      !
      ! Bypass the second call to grid_ahfind()
      ! and compute radius and center in same call:
      !
      logical      BYPASSSECOND
      parameter (  BYPASSSECOND = .true. )
      !
      ! Prescribe motion of hole centers:
      !
      logical      HOLECENTERS
      parameter (  HOLECENTERS  = .false. )
      !
      ! Prevent excision region from shrinking too much
      !
      logical      LIMITSHRINKAGE
      parameter (  LIMITSHRINKAGE = .true. )
      real(kind=8) MIN_AH_RADIUS 
      !
      ! Try to merge distinct regions if "close enough"?
      !     (set threshold 0< <=1.0)
      !
      logical      forceclosemerge
      parameter (  forceclosemerge = .true. )
      real(kind=8) overlap, BHdistance, oversing
      real(kind=8) MERGETHRESH
      parameter (  MERGETHRESH     = 0.0000000001d0 )
      !
      ! Merge once regions are less than this distance from each other:
      !
      real(kind=8) MERGETHRESHDistance
c      parameter (  MERGETHRESHDistance = 3.0d0 )

      !
      ! For debugging purposes,
      ! Force hole #2 to be identical to hole #1:
      !    (modulo reflection to negative x)
      !
      logical      forcesymmetric
      parameter (  forcesymmetric = .false. )

      logical      ltrace
      parameter (  ltrace  = .false. )
      logical      ltrace2
      parameter (  ltrace2 = .false. )
      logical      ltrace3
      parameter (  ltrace3 = .false. )
      ! Field number to be output if tracing is on:
      integer     outfnum,           gridnum
      parameter ( outfnum  = 1 )

            
      gridnum=19
      if(ltrace.and.grid_is_local(gridnum)) then
            call field_out3d(q(gfunc_pointer(outfnum,gridnum)),
     *             gr_t(gridnum),'BHpre',
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             gr_nx(gridnum),gr_ny(gridnum),gr_nz(gridnum),myid)
            call field_out3d(q(gfunc_pointer(outfnum+1,gridnum)),
     *             gr_t(gridnum),'BHpreADV',
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             gr_nx(gridnum),gr_ny(gridnum),gr_nz(gridnum),myid)
      end if
      
      ! set the merge distance as a function of the masses
      MERGETHRESHDistance = 1.3*(bh1_mass + bh2_mass)
c      print*, "MERGETHRESHDistance", MERGETHRESHDistance

      ! set the minimum apparent horizon radius if limiting
      MIN_AH_RADIUS  = 0.7*(bh1_mass + bh2_mass)

      ! If not using the mask, just return
      if (use_mask .le. 0) return

      gi               = level_return_start(level)
      time             = grid_return_time( gi )
      hi               = grid_return_h(    gi )

      !parameter to avoid problems with tiny non-overlaps    
      !small_delta = hxyz0/2**level
      small_delta = hi

      !
      ! Should we bother to mask? Yes if:
      !     (1) the finest level has changed, then remask
      !         (this includes at the initial time)
      !     (2) we are using dynamical excision and
      !         it has been a certain number of iterations
      !

      !write(*,*) 'mask_leve/level: ',mask_level, level
      !write(*,*) 'mask_usemin,mask_minfield: ',mask_usemin,mask_minfield
      !write(*,*) 'mask_period: ', mask_period

      levelchanged   = mask_level .ne. level
      mask_iteration = mask_iteration + 1 

      if (.not. levelchanged) then
         if (ltrace .and. myid.eq.master)
     *   write(*,99)myid,'level_bhmask_local: Level not changed: ',level
         if ( abs(findhorizon) .eq. 0             .and.
     *        abs(bh1_velx) .lt. SMALLNUMBER .and.
     *        abs(bh2_velx) .lt. SMALLNUMBER      ) then
            if (ltrace .and. myid.eq.master)
     *      write(*,99)myid,'level_bhmask_local: Static mask, returning'
            return
         else
            ! If we have not progress a certain number of iterations,
            ! then return
            if( mod(mask_iteration,mask_period).ne.0 ) then
               if (ltrace .and. myid.eq.master)
     *         write(*,99)myid,'level_bhmask_local:Notye',mask_iteration
               return
            else
               if (ltrace .and. myid.eq.master)
     *         write(*,99)myid,'level_bhmask_local: Time to remask...'
            end if
         end if
      else
         if (ltrace .and. myid.eq.master)
     *   write(*,99)myid,'level_bhmask_local: Level changed: ',level
         if (level .gt. mask_level) then
            !
            ! If level has increased then the coords of the mask have changed
            ! because the new level also masks the fine grid point
            ! between the coordinate of the mask and the next coarse gridpoint.
            !
            if (ltrace .and. myid.eq.master)
     *      write(*,99)myid,'level_bhmask_local: Level increased:',level
            do kk = 1, max_num_masks
               if (bh_true(kk)) then
                  do i = 1, 6, 2
                     mask_coords(i,  kk)  = mask_coords(i,  kk)-hi
                     mask_coords(i+1,kk)  = mask_coords(i+1,kk)+hi
                  end do
               end if
              !
              if (ltrace .and. myid.eq.master) then
                 write(*,99)myid,'level_bhmask_local:AChanged Hole#:',kk
                 do i = 1, 6, 2
                    write(*,98)myid,'       mask_coords(i,kk)=',
     *                             mask_coords(i,kk),mask_coords(i+1,kk)
                 end do
              end if
              !
            end do
         end if
         !
         ! Update the new finest level on which we are masking:
         mask_level = level
      end if

      if (ltrace .and. myid.eq.master)
     *        write(*,98)myid,'level_bhmask_local: Masking........',time

      if (ltrace2) then
         write(*,99)myid,'level_bhmask_local: (Re)Locating Horizon(s)'
         write(*,99)myid,'level_bhmask_local:       level= ',level
         write(*,98)myid,'level_bhmask_local:       time = ',time
         write(*,99)myid,'level_bhmask_local:mask_period = ',mask_period
         write(*,99)myid,'level_bhmask_local: mask_itera = ',
     *                mask_iteration
      end if

      !
      ! Store previous mask values so that we can:
      !    --output the mask only when it has changed
      !    --restrict it to moving at most one grid point
      !
      do kk = 1, max_num_masks
         do i = 1, 3
            old_mask_center(i,kk) = mask_center(i,kk)
         end do
         old_mask_radius(1,kk)    = mask_radius(1,kk)
         ! 
         !  Store using new, rectangular style:
         ! 
         do i = 1, 6
            old_mask_bbox(i,kk) = mask_bbox(i,kk)
         end do
         !
         ! If no hole existed, then we do not need to restrict motion:
         if (bh_true(kk)) then 
            if (ltrace .and. myid.eq.master)
     *        write(*,99)myid,'level_bhmask_local:    Previous hole',kk
            noprevioushole(kk)    = .false.
         else
            if (ltrace .and. myid.eq.master)
     *        write(*,99)myid,'level_bhmask_local: No previous hole',kk
            noprevioushole(kk)    = .true.
         end if
      end do

      !
      ! Initialize guess if no hole has been found:
      !   (if hole(s) have been found, use that for guess)
      !
      if (mask_usemin .eq. 0) then
         !
         ! Use user input:
         !
         if (.not. bh_true(1)) then
            if (ltrace2)write(*,99) myid, ' Initializing guess hole 1'
            mask_center(1,1) = bh1_x0
            mask_center(2,1) = bh1_y0
            mask_center(3,1) = bh1_z0
            mask_radius(1,1) = bh1_exc_rad
            !
            mask_bbox(1,1)   = bh1_x0 - bh1_exc_rad
            mask_bbox(2,1)   = bh1_x0 + bh1_exc_rad
            mask_bbox(3,1)   = bh1_y0 - bh1_exc_rad
            mask_bbox(4,1)   = bh1_y0 + bh1_exc_rad
            mask_bbox(5,1)   = bh1_z0 - bh1_exc_rad
            mask_bbox(6,1)   = bh1_z0 + bh1_exc_rad
         end if
         if (.not. bh_true(2)) then
            if (ltrace2)write(*,99) myid, ' Initializing guess hole 2'
            mask_center(1,2) = bh2_x0
            mask_center(2,2) = bh2_y0
            mask_center(3,2) = bh2_z0
            mask_radius(1,2) = bh2_exc_rad
            !
            mask_bbox(1,2)   = bh2_x0 - bh1_exc_rad
            mask_bbox(2,2)   = bh2_x0 + bh1_exc_rad
            mask_bbox(3,2)   = bh2_y0 - bh1_exc_rad
            mask_bbox(4,2)   = bh2_y0 + bh1_exc_rad
            mask_bbox(5,2)   = bh2_z0 - bh1_exc_rad
            mask_bbox(6,2)   = bh2_z0 + bh1_exc_rad
         end if
         !
         if (abs(findhorizon).eq.2) then
            k = 1
            if (.not. bh_true(k)) then
               if(ltrace)write(*,99)myid,'level_bhmask: Initting AH',k
               !
               horizon_center(1,k) = bh1_x0
               horizon_center(2,k) = bh1_y0
               horizon_center(3,k) = bh1_z0
               call load_scal1D(q(horizon_R(k)),bh1_exc_rad,
     *                           horizon_ntheta*horizon_nphi)
            else
               if(ltrace)write(*,96)myid,'level_bhmask: Growing AH',
     *                     k,horizon_growth
               ! If hole already found, allow for possibility that
               ! hole grows, so start with something marginally bigger:
               call vec_scalmult1D(q(horizon_R(k)),horizon_growth,
     *                           horizon_ntheta*horizon_nphi)
            end if
            k = 2
            if (.not. bh_true(k)) then
               if(ltrace)write(*,99)myid,'level_bhmask: Initting AH',k
               !
               horizon_center(1,k) = bh2_x0
               horizon_center(2,k) = bh2_y0
               horizon_center(3,k) = bh2_z0
               call load_scal1D(q(horizon_R(k)),bh2_exc_rad,
     *                           horizon_ntheta*horizon_nphi)
            else
               if(ltrace)write(*,96)myid,'level_bhmask: Growing AH',
     *                     k,horizon_growth
               ! If hole already found, allow for possibility that
               ! hole grows, so start with something marginally bigger:
               call vec_scalmult1D(q(horizon_R(k)),horizon_growth,
     *                           horizon_ntheta*horizon_nphi)
            end if
         end if
         !
      else if (mask_usemin .eq. 1) then
         !
         ! Use minimum of a field:
         !
         if (ltrace)
     *   write(*,*)'level_bhmask: Using minimum of a field for guess',
     *         mask_minfield
         do kk = 1, max_num_masks
            if (.not. bh_true(kk)) then
               !
               ! Find minimum in some field (say lapse):
               !    (1) Find local min on finest level
               !    (2) Find global min across all procs
               !    (3) Broadcast this info from master processor
               !
               loc_min(kk) = LARGENUMBER
               gi = level_return_start(level)
  15           if (grid_return_existence(gi)) then
                  if (grid_is_local(gi)) then
                     if(ltrace2)write(*,99)myid,'call grid_minfind: ',gi
                     call grid_minfind(gi, mask_minfield, loc_min(kk),
     *                                     tmppoint(1,kk),kk           )
                   end if
                   gi = grid_return_sibling(gi)
                   goto 15
               end if
               if (ltrace) then
               write(*,99)myid,'level_bhmask: kk:         ',kk
               write(*,98)myid,'level_bhmask: loc_min(kk):',loc_min(kk)
               write(*,98)myid,'level_bhmask: pnt(1,kk):',tmppoint(1,kk)
               write(*,98)myid,'level_bhmask: pnt(2,kk):',tmppoint(2,kk)
               write(*,98)myid,'level_bhmask: pnt(3,kk):',tmppoint(3,kk)
               end if
               !
               call MPI_Reduce(loc_min(kk),glob_min(kk),1,
     *          MPI_DOUBLE_PRECISION,MPI_MIN,master,MPI_COMM_WORLD,ierr)
               call MPI_BCAST(glob_min(kk), 1, MPI_DOUBLE_PRECISION,
     &                                  master, MPI_COMM_WORLD, ierr)
               ! All procs have the global min, now figure out which
               ! proc has it so we can get coordinates:
               if (double_equal(loc_min(kk),glob_min(kk))) then
                  call MPI_Reduce(myid,proc_min(kk),1,MPI_INTEGER,
     *                         MPI_MIN,master,MPI_COMM_WORLD,ierr)
               else
                  call MPI_Reduce(numprocs,proc_min(kk),1,MPI_INTEGER,
     *                         MPI_MIN,master,MPI_COMM_WORLD,ierr)
               end if
               ! Let everyone know that proc_min will broadcast coords:
               call MPI_BCAST(proc_min(kk), 1, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
               ! Everyone now gets the coords of the min from proc_min:
               call MPI_BCAST(tmppoint(1,kk), 3, MPI_DOUBLE_PRECISION,
     &                            proc_min(kk), MPI_COMM_WORLD, ierr)
               if (ltrace) then
               write(*,99)myid,'level_bhmask: kk:         ',kk
               write(*,99)myid,'level_bhmask: proc_min:   ',proc_min(kk)
               write(*,98)myid,'level_bhmask:glob_min(kk):',glob_min(kk)
               write(*,98)myid,'level_bhmask: pnt(1,kk):',tmppoint(1,kk)
               write(*,98)myid,'level_bhmask: pnt(2,kk):',tmppoint(2,kk)
               write(*,98)myid,'level_bhmask: pnt(3,kk):',tmppoint(3,kk)
               end if
               !
               do i = 1, 3
                  mask_center(i,kk) = tmppoint(i,kk)
               end do
               if (kk.eq.1) then
                  mask_radius(1,kk) = bh1_exc_rad
               else if (kk.eq.2) then
                  mask_radius(1,kk) = bh2_exc_rad
               end if
               !
               mask_bbox(1,1) = tmppoint(1,1) - bh1_exc_rad
               mask_bbox(2,1) = tmppoint(1,1) + bh1_exc_rad
               mask_bbox(3,1) = tmppoint(2,1) - bh1_exc_rad
               mask_bbox(4,1) = tmppoint(2,1) + bh1_exc_rad
               mask_bbox(5,1) = tmppoint(3,1) - bh1_exc_rad
               mask_bbox(6,1) = tmppoint(3,1) + bh1_exc_rad
               !
               mask_bbox(1,2) = tmppoint(1,2) - bh2_exc_rad
               mask_bbox(2,2) = tmppoint(1,2) + bh2_exc_rad
               mask_bbox(3,2) = tmppoint(2,2) - bh2_exc_rad
               mask_bbox(4,2) = tmppoint(2,2) + bh2_exc_rad
               mask_bbox(5,2) = tmppoint(3,2) - bh2_exc_rad
               mask_bbox(6,2) = tmppoint(3,2) + bh2_exc_rad
               !
               if (myid.eq.master) then
                 open(unit=10,STATUS='UNKNOWN')
                 write(10,"(9E15.7)") time,
     *           glob_min(1), tmppoint(1,1),tmppoint(2,1),tmppoint(3,1),
     *           glob_min(2), tmppoint(1,2),tmppoint(2,2),tmppoint(3,2)
               end if
               !
            end if
         end do
      else
         write(*,*)'level_bhmask: Unknown val mask_usemin: ',mask_usemin
      end if

      !
      ! Determine location of masked regions:
      !
      if (abs(findhorizon).eq.1) then
         !........................................................
         !      I N F L O W    F I N D E R  
         !........................................................
         if (ltrace.and.myid.eq.master)
     *   write(*,99) myid, ' Looking for a horizon...'
         !
         ! Iterate on the center:
         !
         iter = 0
 23      if (iter .lt. MAXITERS) then
            iter = iter + 1
            if (ltrace.and.myid.eq.master) write(*,99)
     *         myid, 'level_bhmask:  Iteration: ',iter
            
         !
         ! Use a horizon finder for each black hole
         !
         !
         if (ltrace2)write(*,99) myid, ' Initialize bounds local'
         do kk = 1, max_num_masks
            bh_true(kk)       = .false.
            ! Initialize local bounds:
            local_left(1,kk)  = maxx0
            local_right(1,kk) = minx0
            local_left(2,kk)  = maxy0
            local_right(2,kk) = miny0
            local_left(3,kk)  = maxz0
            local_right(3,kk) = minz0
         end do
         !
         ! the first call to grid_ahfind is to update the centers
         !     (grid_ahfind loops over both regions)
         !     NB: now loops over levels, starting from finest
         !         continuing as long as coarser levels are time
         !         aligned:
         lev_i = level
  19     if (level_time_aligned(lev_i,level) ) then
            if (myid.eq.master)
     *         write(*,99)myid, 'A: Time aligned lev_i:',lev_i
            gi = level_return_start(lev_i)
  20        if (grid_return_existence(gi)) then
              if (grid_is_local(gi)) then
                if (ltrace2)write(*,99)myid, ' calling grid_ahfind: ',gi
                ! could exclude grids which are fully covered by finer level
                call grid_ahfind(gi)
                if (ltrace2)
     *   call field_out3d(q(gfunc_pointer(1,gi)),gr_t(gi),
     *             'bhasmk_exp',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             gr_nx(gi),gr_ny(gi),gr_nz(gi),myid)
              end if
              gi = grid_return_sibling(gi)
              goto 20
            end if
            !
            lev_i = lev_i - 1
            if (lev_i .ge. 0) goto 19
            !
         end if
         !
         ! If we merged holes already, do not allow a second hole:
         !
         if (forcedmerge) bh_true(2) = .false.
         !
         if (ltrace2) then
            do kk = 1, max_num_masks
               write(*,99)myid,' Before bh_true(kk):   ',kk,bh_true(kk)
            end do
         end if
         do kk = 1, max_num_masks
            tmp_bh_true(kk) = bh_true(kk)
         end do
         call MPI_Reduce(tmp_bh_true,bh_true, max_num_masks,
     *         MPI_LOGICAL,MPI_LOR,master,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(bh_true, max_num_masks, MPI_LOGICAL,
     &                                  master, MPI_COMM_WORLD, ierr)
         if (ltrace2) then
            do kk = 1, max_num_masks
               write(*,99)myid,' After bh_true(kk):   ',kk,bh_true(kk)
            end do
         end if
         !
         ! if any inflow boundaries have been found then compute the
         ! center and the radius. 
         !
         do kk = 1, max_num_masks
            if (bh_true(kk)) then
               ! this proc should now have bounds local to this proc
               ! need to use MPI to find global bounds:
               if (ltrace3) then
                   write(*,99)myid,'    Local  bnds for kk = ',kk
                   do i = 1, 3
                     write(*,98)myid,'    local_left/right(i,kk) =',
     *               local_left(i,kk),local_right(i,kk)
                   end do
               end if
               !
               ! compute the global minimum/maximum distance from the
               ! old center to the inflow boundary
               if (ltrace3)write(*,99)myid,'    Calling MPI_Reduce'
               call MPI_Reduce(local_left(1,kk),mask_left(1,kk), 3,
     *          MPI_DOUBLE_PRECISION,MPI_MIN,master,MPI_COMM_WORLD,ierr)
               call MPI_Reduce(local_right(1,kk),mask_right(1,kk),3,
     *          MPI_DOUBLE_PRECISION,MPI_MAX,master,MPI_COMM_WORLD,ierr)
               call MPI_BCAST(mask_left(1,kk),3,
     *                  MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
               call MPI_BCAST(mask_right(1,kk),3,
     *                  MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
               !
               if (ltrace3 .or. (ltrace2.and. myid.eq.master)) then
                   write(*,99)myid,'    Global bnds for kk = ',kk
                   do i = 1, 3
                      write(*,98)myid,'    mask_left/right(i,kk) =',
     *                   mask_left(i,kk),mask_right(i,kk)
                   end do
               end if
               ! check that the center is outside the inflow boundary
               if ( (mask_center(1,kk) .lt. mask_left(1,kk)) .or.
     *              (mask_center(2,kk) .lt. mask_left(2,kk)) .or.
     *              (mask_center(3,kk) .lt. mask_left(3,kk)) .or.
     *              (mask_center(1,kk) .gt. mask_right(1,kk)) .or.
     *              (mask_center(2,kk) .gt. mask_right(2,kk)) .or.
     *              (mask_center(3,kk) .gt. mask_right(3,kk))  ) then
                   ! We need to remove this masked region:
                   bh_true(kk) = .false.
                  if (ltrace2) write(*,99) myid,
     *                '  Removing region bcause center out of bds',kk
                  if (ltrace2) then
                      write(*,98)myid,'1l/c/r:',
     *                mask_left(1,kk),mask_center(1,kk),mask_right(1,kk)
                      write(*,98)myid,'2l/c/r:',
     *                mask_left(2,kk),mask_center(2,kk),mask_right(2,kk)
                      write(*,98)myid,'3l/c/r:',
     *                mask_left(3,kk),mask_center(3,kk),mask_right(3,kk)
                  end if
               else if (ltrace2) then
                  write(*,99)myid,'level_bhmask_local:Cntr w/in bnds',kk
               end if
               !
               ! Compute the updated center
               !
              if (.not. HOLECENTERS) then
              mask_center(1,kk)=0.5d0*(mask_left(1,kk)+mask_right(1,kk))
              mask_center(2,kk)=0.5d0*(mask_left(2,kk)+mask_right(2,kk))
              mask_center(3,kk)=0.5d0*(mask_left(3,kk)+mask_right(3,kk))
              else
                 ! Have centers move along some curve in time:
                 !motion            = 0.1d0 * time + 0.15d0 * time**2
                 motion            = 0.13d0 * time - 0.05d0 * time**2
                 if (kk.eq.2) motion = -motion
                 mask_center(1,kk) = mask_center(1,kk) + motion
                 mask_center(2,kk) = 0.0d0
                 mask_center(3,kk) = 0.0d0
                 if (myid.eq.master.and.kk.eq.1) then 
                   write(*,*)'motion,center,t: ',
     *                         motion,mask_center(1,kk),time
                 end if
              end if
               !
               ! Only the root has the correct answer, so let us bcast:
               if (ltrace3)write(*,99)myid,'    Calling MPI_BCAST'
               call MPI_BCAST(mask_center(1,kk),3,
     *                  MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
               !
               if (ltrace2) then
                    do i = 1, 3
                       write(*,98)myid,'     center: ',mask_center(i,kk)
                    end do
               end if
               !
            else
               if (ltrace2)write(*,99)myid,'    No mask in region: ',kk
            end if
            !
         end do

         !
         ! We have the old centers and the new ones
         ! Check if the difference is bigger than 1 grid point, 
         !    otherwise we can keep the old one
         !
         do kk = 1, max_num_masks
            if (bh_true(kk)) then
              mindistance = 3.01*hi**2
!              mindistance = 1.51*hi**2
              if(forcedmerge) mindistance = 3.001*hi**2

              distance    =
     *                    (mask_center(1,kk) - old_mask_center(1,kk))**2
     *                   +(mask_center(2,kk) - old_mask_center(2,kk))**2
     *                   +(mask_center(3,kk) - old_mask_center(3,kk))**2
              if (distance .lt. mindistance) then
                 if (ltrace2)write(*,99)myid,'    Using old center'
                 do pp=1,3
                  mask_center(pp,kk) = old_mask_center(pp,kk)
                 end do
              end if
              !
              ! Only the root has the correct answer, so let us bcast:
              if (ltrace3)write(*,99)myid,'    Calling MPI_BCAST'
              call MPI_BCAST(mask_center(1,kk),3,
     *                 MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
              !
              if (ltrace2) then
                 do i = 1, 3
                    write(*,98)myid,'     center: ',mask_center(i,kk)
                 end do
              end if
              !
            else
               if (ltrace2)write(*,99)myid,'    No mask in region: ',kk
            end if
            !
         end do
         
         !
            !
            if (ltrace.and.myid.eq.master) then
               do kk = 1, max_num_masks
                  if (bh_true(kk)) then
                     write(*,99)myid,'     iter/kk: ', iter, kk
                     do i = 1, 3
                 write(*,98)myid,'center:',mask_center(i,kk)
                 !write(*,98)myid,'l/r:',mask_left(i,kk),mask_right(i,kk)
                     end do
                  end if
               end do
            end if
            ! go back and try another shot at finding the center
            goto 23
         end if

         !
         ! At this point we have the centers for 
         ! num_masks mask regions, but now we have to determine their radii
         !     (NB: we use tmp_num_masks because num_masks gets changed
         !      by the call to grid_ahfind())
         !
         tmp_num_masks = num_masks
         do kk = 1, max_num_masks
            tmp_bh_true(kk) = bh_true(kk)
         end do
         !
         ! Skip second call to grid_ahfind()?
         !
         if (BYPASSSECOND) goto 33
         !if (.true.) goto 33
         !
         !
         if(ltrace2)write(*,99)myid,'Finding radii for regions '
         do kk = 1, max_num_masks
               ! Initialize local bounds:
               local_left(1,kk)  = maxx0
               local_right(1,kk) = minx0
               local_left(2,kk)  = maxy0
               local_right(2,kk) = miny0
               local_left(3,kk)  = maxz0
               local_right(3,kk) = minz0
         end do
         !
         lev_i = level
  29     if (level_time_aligned(lev_i,level) ) then
            if (myid.eq.master)
     *         write(*,99)myid, 'B: Time aligned lev_i:',lev_i
            gi = level_return_start(lev_i)
  30        if (grid_return_existence(gi)) then
               if (grid_is_local(gi)) then
                  if(ltrace2)write(*,99)myid,'callng grid_ahfind:',gi
                  call grid_ahfind(gi)
                  if (ltrace2)
     *   call field_out3d(q(gfunc_pointer(1,gi)),gr_t(gi),
     *             'bhasmk_expB',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             gr_nx(gi),gr_ny(gi),gr_nz(gi),myid)
               end if
               gi = grid_return_sibling(gi)
               goto 30
            end if
            !
            lev_i = lev_i - 1
            if (lev_i .ge. 0) goto 29
            !
         end if
         ! this proc should now have bounds local to this proc
         ! need to use MPI to find global bounds:
         if (ltrace3) then
            do kk= 1, max_num_masks
               if (tmp_bh_true(kk)) then
                    do i = 1, 3
                       write(*,98)myid,'    local_left/right(i,kk) =',
     *                             local_left(i,kk),local_right(i,kk)
                    end do
               end if
            end do
         end if
         ! compute the global minimum/maximum distance from the
         ! old center to the inflow boundary
         if (ltrace3)write(*,99)myid,'    Calling MPI_Reduce'
         call MPI_Reduce(local_left(1,1),mask_left(1,1),3*max_num_masks,
     *          MPI_DOUBLE_PRECISION,MPI_MIN,master,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(mask_left(1,1), 3*max_num_masks, 
     *                  MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
         call MPI_Reduce(local_right(1,1),mask_right(1,1),
     *               3*max_num_masks,
     *          MPI_DOUBLE_PRECISION,MPI_MAX,master,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(mask_right(1,1), 3*max_num_masks, 
     *                  MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
         !
 33      if (ltrace3 ) then
         !if (ltrace3 .and. myid.eq.master) then
            do kk= 1, max_num_masks
               if (tmp_bh_true(kk)) then
                   write(*,99)myid,'    Global Bounds for Region #: ',kk
                   do i = 1, 3
                      write(*,98)myid,'    mask_left/right(i,kk) =',
     *                mask_left(i,kk),mask_right(i,kk)
                   end do
               end if
            end do
         end if
         !
         do kk = 1, max_num_masks
            if (tmp_bh_true(kk)) then
                if(kk.eq.1) weight = bh1_mass
                if(kk.eq.2) weight = bh2_mass
                if(forcedmerge) weight = bh1_mass + bh2_mass
               !
               ! Use new rectangular style:
               !    (want to excise within a shrunken region
               !      w/r/t mask_left/right)
               !
               do i = 1, 6, 2
                  ii = (i+1)/2
               span = 0.5d0*(mask_right(ii,kk) - mask_left(ii,kk))
               if(span.lt.weight*1.4.and.local_time.gt.2)
     *            span = .7*weight/SHRINK

                  mask_bbox(i,  kk) =
     *                     0.5d0*( mask_right(ii,kk) + mask_left(ii,kk))
     *          - SHRINK * span
                  mask_bbox(i+1,kk) =
     *                     0.5d0*( mask_right(ii,kk) + mask_left(ii,kk))
     *          + SHRINK * span


c                  mask_bbox(i,  kk) =
c     *                     0.5d0*( mask_right(ii,kk) + mask_left(ii,kk))
c     *          - SHRINK * 0.5d0*( mask_right(ii,kk) - mask_left(ii,kk))
c                  mask_bbox(i+1,kk) =
c     *                     0.5d0*( mask_right(ii,kk) + mask_left(ii,kk))
c     *          + SHRINK * 0.5d0*( mask_right(ii,kk) - mask_left(ii,kk))
c                  if (ltrace2 .and. myid.eq.master)
c                  !if (ltrace .and. myid.eq.master)
c     *                 write(*,*) 'mask_bbox:',i,
c     *                    mask_bbox(i,kk), mask_bbox(i+1,kk)
               end do
               !
               ! The center is the same as before 
               ! and now we compute the radius         
               rad_left_x  = abs(mask_left(1,kk)  - mask_center(1,kk))
               rad_right_x = abs(mask_right(1,kk) - mask_center(1,kk))
               rad_left_y  = abs(mask_left(2,kk)  - mask_center(2,kk))
               rad_right_y = abs(mask_right(2,kk) - mask_center(2,kk))
               rad_left_z  = abs(mask_left(3,kk)  - mask_center(3,kk))
               rad_right_z = abs(mask_right(3,kk) - mask_center(3,kk))
               !
               mask_radius(1,kk) = min(rad_left_x, rad_right_x)
               mask_radius(2,kk) = min(rad_left_y, rad_right_y)
               mask_radius(3,kk) = min(rad_left_z, rad_right_z)
               !
               if (ltrace3) then
                    do i = 1, 3
                        write(*,98)myid,' radius in each direction =',
     *              mask_radius(i,kk)
                     end do
               end if
               !mask_radius(1,kk) = 0.75d0*min(mask_radius(1,kk), 
               mask_radius(1,kk) = 0.5d0*min(mask_radius(1,kk), 
     &                            mask_radius(2,kk),mask_radius(3,kk))
               call MPI_BCAST(mask_radius(1,kk), 1, 
     *                  MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
               if(ltrace2)write(*,98)myid,'glb radius',mask_radius(1,kk)
               !
            else
               !
               if(ltrace2)write(*,99)myid,'No masking   for region: ',kk
               !
            end if
         end do  
         !
         ! Restore num_masks to proper value:
         num_masks = tmp_num_masks   
         do kk = 1, max_num_masks
            bh_true(kk) = tmp_bh_true(kk)
         end do
      else if (abs(findhorizon).eq.0) then
         !........................................................
         !      S T A T I C     M A S K
         !........................................................
         if (.not.forcedmerge) then
         !
         ! Move masks linearly:
         !
         num_masks        = nbholes
         do kk = 1, nbholes
            bh_true(kk) = .true.
         end do
         !
         motion           = bh1_velx * time
         mask_center(1,1) = bh1_x0 + motion
         mask_center(2,1) = bh1_y0
         mask_center(3,1) = bh1_z0
         mask_radius(1,1) = bh1_exc_rad
         !
         mask_bbox(1,1)   = bh1_x0 - bh1_exc_rad + motion
         mask_bbox(2,1)   = bh1_x0 + bh1_exc_rad + motion
         mask_bbox(3,1)   = bh1_y0 - bh1_exc_rad
         mask_bbox(4,1)   = bh1_y0 + bh1_exc_rad
         mask_bbox(5,1)   = bh1_z0 - bh1_exc_rad
         mask_bbox(6,1)   = bh1_z0 + bh1_exc_rad
         !
         !
         if (ltrace) write(*,98) myid, '   center1: ',mask_center(1,1)
      !write(*,60) ' level_bhmask_local: t/h/motion: ',time, hi, motion
 60   format(A,3F10.5)
         !
         motion           = bh2_velx * time
         mask_center(1,2) = bh2_x0 + motion
         mask_center(2,2) = bh2_y0
         mask_center(3,2) = bh2_z0
         mask_radius(1,2) = bh2_exc_rad
         !
         mask_bbox(1,2)   = bh2_x0 - bh2_exc_rad + motion
         mask_bbox(2,2)   = bh2_x0 + bh2_exc_rad + motion
         mask_bbox(3,2)   = bh2_y0 - bh2_exc_rad
         mask_bbox(4,2)   = bh2_y0 + bh2_exc_rad
         mask_bbox(5,2)   = bh2_z0 - bh2_exc_rad
         mask_bbox(6,2)   = bh2_z0 + bh2_exc_rad
         !
         !
         if (.true.) then
            ! Circular motion:
            !
            omega            = 3.0d0
            !
            motionx          = -bh1_velx * sin( omega*time+3.14d0)
     *                         * abs(bh1_x0) * time
            motiony          = +bh1_velx * cos( omega*time+3.14d0)
     *                         * abs(bh1_x0) * time
            !write(*,*) '1: motionx/y: ',motionx,motiony
            !
            !write(*,*) '1: bbox1/2: ',mask_bbox(1,1),mask_bbox(2,1)
            mask_bbox(1,1)   = bh1_x0 - bh1_exc_rad + motionx
            mask_bbox(2,1)   = bh1_x0 + bh1_exc_rad + motionx
            mask_bbox(3,1)   = bh1_y0 - bh1_exc_rad + motiony
            mask_bbox(4,1)   = bh1_y0 + bh1_exc_rad + motiony
            mask_bbox(5,1)   = bh1_z0 - bh1_exc_rad
            mask_bbox(6,1)   = bh1_z0 + bh1_exc_rad
            !
            mask_center(1,1) = 0.5d0*(mask_bbox(1,1)+mask_bbox(2,1))
            mask_center(2,1) = 0.5d0*(mask_bbox(3,1)+mask_bbox(4,1))
            mask_center(3,1) = bh1_z0
            mask_radius(1,1) = bh1_exc_rad
            !
            !
c           motionx          = -bh2_velx * sin( omega*time)
c    *                         * abs(bh2_x0) * time
c           motiony          = +bh2_velx * cos( omega*time)
c    *                         * abs(bh2_x0) * time
            ! Ensure holes move diametrically opposite on circle:
            motionx          = -motionx
            motiony          = -motiony
            !
            mask_bbox(1,2)   = bh2_x0 - bh2_exc_rad + motionx
            mask_bbox(2,2)   = bh2_x0 + bh2_exc_rad + motionx
            mask_bbox(3,2)   = bh2_y0 - bh2_exc_rad + motiony
            mask_bbox(4,2)   = bh2_y0 + bh2_exc_rad + motiony
            mask_bbox(5,2)   = bh2_z0 - bh2_exc_rad
            mask_bbox(6,2)   = bh2_z0 + bh2_exc_rad
            !
            mask_center(1,2) = 0.5d0*(mask_bbox(1,2)+mask_bbox(2,2))
            mask_center(2,2) = 0.5d0*(mask_bbox(3,2)+mask_bbox(4,2))
            mask_center(3,2) = bh2_z0
            mask_radius(1,2) = bh2_exc_rad
            !
            !write(*,*) '1: bbox1/2: ',mask_bbox(1,1),mask_bbox(2,1)
         end if
         !
         if (ltrace) write(*,98) myid, '   center2: ',mask_center(1,2)
         end if
         !
      else if (abs(findhorizon).eq.2) then
         !........................................................
         !      H O R I Z O N     F I N D E R
         !........................................................
         if (myid.eq.master.and.ltrace)write(*,*)'Using Horizon finder'
         call horizon_evolve()
         !
      else if (abs(findhorizon).eq.3) then
         !........................................................
         !      M I J A N ' S     F I N D E R
         !........................................................
          if (myid.eq.master.and.ltrace)write(*,*)'Using Mijans finder'
            call mijans_horizon(level)
      else if (abs(findhorizon).eq.4) then
      ! This is just a fixed excision size for testing
        bh_true(1) = .true.
        mask_bbox(1,1) = -1.0d0
        mask_bbox(3,1) = -1.0d0
        mask_bbox(5,1) = -1.0d0
        mask_bbox(2,1) = 1.0d0
        mask_bbox(4,1) = 1.0d0
        mask_bbox(6,1) = 1.0d0
        mask_center(1,1) = 0.0d0
        mask_center(2,1) = 0.0d0
        mask_center(3,1) = 0.0d0
      else if (abs(findhorizon).eq.5) then
         !........................................................
         !      I N F L O W    F I N D E R  
         !........................................................
         if (ltrace.and.myid.eq.master)
     *   write(*,99) myid, ' Looking for a horizon...'
         !
         ! Iterate on the center:
         !
         iter = 0
 123      if (iter .lt. MAXITERS) then
            iter = iter + 1
            if (ltrace.and.myid.eq.master) write(*,99)
     *         myid, 'level_bhmask:  Iteration: ',iter
            
         !
         ! Use a horizon finder for each black hole
         !
         !
         if (ltrace2)write(*,99) myid, ' Initialize bounds local'
         do kk = 1, max_num_masks
            bh_true(kk)       = .false.
            ! Initialize local bounds:
            local_left(1,kk)  = maxx0
            local_right(1,kk) = minx0
            local_left(2,kk)  = maxy0
            local_right(2,kk) = miny0
            local_left(3,kk)  = maxz0
            local_right(3,kk) = minz0
         end do
         !
         ! the first call to grid_ahfind is to update the centers
         !     (grid_ahfind loops over both regions)
         !     NB: now loops over levels, starting from finest
         !         continuing as long as coarser levels are time
         !         aligned:
         lev_i = level
  119     if (level_time_aligned(lev_i,level) ) then
            if (myid.eq.master)
     *         write(*,99)myid, 'A: Time aligned lev_i:',lev_i
            gi = level_return_start(lev_i)
  120        if (grid_return_existence(gi)) then
              if (grid_is_local(gi)) then
                if (ltrace2)write(*,99)myid, ' calling grid_ahfind: ',gi
                ! could exclude grids which are fully covered by finer level
                call grid_ahfind(gi)
                if (ltrace2)
     *   call field_out3d(q(gfunc_pointer(1,gi)),gr_t(gi),
     *             'bhasmk_exp',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             gr_nx(gi),gr_ny(gi),gr_nz(gi),myid)
              end if
              gi = grid_return_sibling(gi)
              goto 120
            end if
            !
            lev_i = lev_i - 1
            if (lev_i .ge. 0) goto 119
            !
         end if
         !
         ! If we merged holes already, do not allow a second hole:
         !
         if (forcedmerge) bh_true(2) = .false.
         !
         if (ltrace2) then
            do kk = 1, max_num_masks
               write(*,99)myid,' Before bh_true(kk):   ',kk,bh_true(kk)
            end do
         end if
         do kk = 1, max_num_masks
            tmp_bh_true(kk) = bh_true(kk)
         end do
         call MPI_Reduce(tmp_bh_true,bh_true, max_num_masks,
     *         MPI_LOGICAL,MPI_LOR,master,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(bh_true, max_num_masks, MPI_LOGICAL,
     &                                  master, MPI_COMM_WORLD, ierr)
         if (ltrace2) then
            do kk = 1, max_num_masks
               write(*,99)myid,' After bh_true(kk):   ',kk,bh_true(kk)
            end do
         end if
         !
         ! if any inflow boundaries have been found then compute the
         ! center and the radius. 
         !
         do kk = 1, max_num_masks
            if (bh_true(kk)) then
               ! this proc should now have bounds local to this proc
               ! need to use MPI to find global bounds:
               if (ltrace3) then
                   write(*,99)myid,'    Local  bnds for kk = ',kk
                   do i = 1, 3
                     write(*,98)myid,'    local_left/right(i,kk) =',
     *               local_left(i,kk),local_right(i,kk)
                   end do
               end if
               !
               ! compute the global minimum/maximum distance from the
               ! old center to the inflow boundary
               if (ltrace3)write(*,99)myid,'    Calling MPI_Reduce'
               call MPI_Reduce(local_left(1,kk),mask_left(1,kk), 3,
     *          MPI_DOUBLE_PRECISION,MPI_MIN,master,MPI_COMM_WORLD,ierr)
               call MPI_Reduce(local_right(1,kk),mask_right(1,kk),3,
     *          MPI_DOUBLE_PRECISION,MPI_MAX,master,MPI_COMM_WORLD,ierr)
               call MPI_BCAST(mask_left(1,kk),3,
     *                  MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
               call MPI_BCAST(mask_right(1,kk),3,
     *                  MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
               !
               if (ltrace3 .or. (ltrace2.and. myid.eq.master)) then
                   write(*,99)myid,'    Global bnds for kk = ',kk
                   do i = 1, 3
                      write(*,98)myid,'    mask_left/right(i,kk) =',
     *                   mask_left(i,kk),mask_right(i,kk)
                   end do
               end if
               ! check that the center is outside the inflow boundary
               if ( (mask_center(1,kk) .lt. mask_left(1,kk)) .or.
     *              (mask_center(2,kk) .lt. mask_left(2,kk)) .or.
     *              (mask_center(3,kk) .lt. mask_left(3,kk)) .or.
     *              (mask_center(1,kk) .gt. mask_right(1,kk)) .or.
     *              (mask_center(2,kk) .gt. mask_right(2,kk)) .or.
     *              (mask_center(3,kk) .gt. mask_right(3,kk))  ) then
                   ! We need to remove this masked region:
                   bh_true(kk) = .false.
                  if (ltrace2) write(*,99) myid,
     *                '  Removing region bcause center out of bds',kk
                  if (ltrace2) then
                      write(*,98)myid,'1l/c/r:',
     *                mask_left(1,kk),mask_center(1,kk),mask_right(1,kk)
                      write(*,98)myid,'2l/c/r:',
     *                mask_left(2,kk),mask_center(2,kk),mask_right(2,kk)
                      write(*,98)myid,'3l/c/r:',
     *                mask_left(3,kk),mask_center(3,kk),mask_right(3,kk)
                  end if
               else if (ltrace2) then
                  write(*,99)myid,'level_bhmask_local:Cntr w/in bnds',kk
               end if
               !
               ! Compute the updated center
               !
              if (.not. HOLECENTERS) then
              mask_center(1,kk)=0.5d0*(mask_left(1,kk)+mask_right(1,kk))
              mask_center(2,kk)=0.5d0*(mask_left(2,kk)+mask_right(2,kk))
              mask_center(3,kk)=0.5d0*(mask_left(3,kk)+mask_right(3,kk))
              else
                 ! Have centers move along some curve in time:
                 !motion            = 0.1d0 * time + 0.15d0 * time**2
                 motion            = 0.13d0 * time - 0.05d0 * time**2
                 if (kk.eq.2) motion = -motion
                 mask_center(1,kk) = mask_center(1,kk) + motion
                 mask_center(2,kk) = 0.0d0
                 mask_center(3,kk) = 0.0d0
                 if (myid.eq.master.and.kk.eq.1) then 
                   write(*,*)'motion,center,t: ',
     *                         motion,mask_center(1,kk),time
                 end if
              end if
               !
               ! Only the root has the correct answer, so let us bcast:
               if (ltrace3)write(*,99)myid,'    Calling MPI_BCAST'
               call MPI_BCAST(mask_center(1,kk),3,
     *                  MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
               !
               if (ltrace2) then
                    do i = 1, 3
                       write(*,98)myid,'     center: ',mask_center(i,kk)
                    end do
               end if
               !
            else
               if (ltrace2)write(*,99)myid,'    No mask in region: ',kk
            end if
            !
         end do

         !
         ! We have the old centers and the new ones
         ! Check if the difference is bigger than 1 grid point, 
         !    otherwise we can keep the old one
         !
         do kk = 1, max_num_masks
            if (bh_true(kk)) then
!              mindistance = 3.01*hi**2
              mindistance = 1.51*hi**2
              if(forcedmerge) mindistance = 3.001*hi**2

              distance    =
     *                    (mask_center(1,kk) - old_mask_center(1,kk))**2
     *                   +(mask_center(2,kk) - old_mask_center(2,kk))**2
     *                   +(mask_center(3,kk) - old_mask_center(3,kk))**2
              if (distance .lt. mindistance) then
                 if (ltrace2)write(*,99)myid,'    Using old center'
                 do pp=1,3
                  mask_center(pp,kk) = old_mask_center(pp,kk)
                 end do
              end if
              !
              ! Only the root has the correct answer, so let us bcast:
              if (ltrace3)write(*,99)myid,'    Calling MPI_BCAST'
              call MPI_BCAST(mask_center(1,kk),3,
     *                 MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
              !
              if (ltrace2) then
                 do i = 1, 3
                    write(*,98)myid,'     center: ',mask_center(i,kk)
                 end do
              end if
              !
            else
               if (ltrace2)write(*,99)myid,'    No mask in region: ',kk
            end if
            !
         end do
         
         !
            !
            if (ltrace.and.myid.eq.master) then
               do kk = 1, max_num_masks
                  if (bh_true(kk)) then
                     write(*,99)myid,'     iter/kk: ', iter, kk
                     do i = 1, 3
                 write(*,98)myid,'center:',mask_center(i,kk)
                 !write(*,98)myid,'l/r:',mask_left(i,kk),mask_right(i,kk)
                     end do
                  end if
               end do
            end if
            ! go back and try another shot at finding the center
            goto 123
         end if

         !
         ! At this point we have the centers for 
         ! num_masks mask regions, but now we have to determine their radii
         !     (NB: we use tmp_num_masks because num_masks gets changed
         !      by the call to grid_ahfind())
         !
         tmp_num_masks = num_masks
         do kk = 1, max_num_masks
            tmp_bh_true(kk) = bh_true(kk)
         end do
         !
         ! Skip second call to grid_ahfind()?
         !
         if (BYPASSSECOND) goto 133
         !if (.true.) goto 33
         !
         !
         if(ltrace2)write(*,99)myid,'Finding radii for regions '
         do kk = 1, max_num_masks
               ! Initialize local bounds:
               local_left(1,kk)  = maxx0
               local_right(1,kk) = minx0
               local_left(2,kk)  = maxy0
               local_right(2,kk) = miny0
               local_left(3,kk)  = maxz0
               local_right(3,kk) = minz0
         end do
         !
         lev_i = level
  129     if (level_time_aligned(lev_i,level) ) then
            if (myid.eq.master)
     *         write(*,99)myid, 'B: Time aligned lev_i:',lev_i
            gi = level_return_start(lev_i)
  130        if (grid_return_existence(gi)) then
               if (grid_is_local(gi)) then
                  if(ltrace2)write(*,99)myid,'callng grid_ahfind:',gi
                  call grid_ahfind(gi)
                  if (ltrace2)
     *   call field_out3d(q(gfunc_pointer(1,gi)),gr_t(gi),
     *             'bhasmk_expB',
     *             gr_minx(gi),gr_maxx(gi),
     *             gr_miny(gi),gr_maxy(gi),
     *             gr_minz(gi),gr_maxz(gi),
     *             gr_nx(gi),gr_ny(gi),gr_nz(gi),myid)
               end if
               gi = grid_return_sibling(gi)
               goto 130
            end if
            !
            lev_i = lev_i - 1
            if (lev_i .ge. 0) goto 129
            !
         end if
         ! this proc should now have bounds local to this proc
         ! need to use MPI to find global bounds:
         if (ltrace3) then
            do kk= 1, max_num_masks
               if (tmp_bh_true(kk)) then
                    do i = 1, 3
                       write(*,98)myid,'    local_left/right(i,kk) =',
     *                             local_left(i,kk),local_right(i,kk)
                    end do
               end if
            end do
         end if
         ! compute the global minimum/maximum distance from the
         ! old center to the inflow boundary
         if (ltrace3)write(*,99)myid,'    Calling MPI_Reduce'
         call MPI_Reduce(local_left(1,1),mask_left(1,1),3*max_num_masks,
     *          MPI_DOUBLE_PRECISION,MPI_MIN,master,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(mask_left(1,1), 3*max_num_masks, 
     *                  MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
         call MPI_Reduce(local_right(1,1),mask_right(1,1),
     *               3*max_num_masks,
     *          MPI_DOUBLE_PRECISION,MPI_MAX,master,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(mask_right(1,1), 3*max_num_masks, 
     *                  MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
         !
 133      if (ltrace3 ) then
         !if (ltrace3 .and. myid.eq.master) then
            do kk= 1, max_num_masks
               if (tmp_bh_true(kk)) then
                   write(*,99)myid,'    Global Bounds for Region #: ',kk
                   do i = 1, 3
                      write(*,98)myid,'    mask_left/right(i,kk) =',
     *                mask_left(i,kk),mask_right(i,kk)
                   end do
               end if
            end do
         end if
         !
         do kk = 1, max_num_masks
            if (tmp_bh_true(kk)) then
               !
               ! Use new rectangular style:
               !    (want to excise within a shrunken region
               !      w/r/t mask_left/right)
               !
               do i = 1, 6, 2
                  ii = (i+1)/2
                  mask_bbox(i,  kk) =
     *                     0.5d0*( mask_right(ii,kk) + mask_left(ii,kk))
     *          - SHRINK * 0.5d0*( mask_right(ii,kk) - mask_left(ii,kk))
                  mask_bbox(i+1,kk) =
     *                     0.5d0*( mask_right(ii,kk) + mask_left(ii,kk))
     *          + SHRINK * 0.5d0*( mask_right(ii,kk) - mask_left(ii,kk))
                  if (ltrace2 .and. myid.eq.master)
                  !if (ltrace .and. myid.eq.master)
     *                 write(*,*) 'mask_bbox:',i,
     *                    mask_bbox(i,kk), mask_bbox(i+1,kk)
               end do
               !
               ! The center is the same as before 
               ! and now we compute the radius         
               rad_left_x  = abs(mask_left(1,kk)  - mask_center(1,kk))
               rad_right_x = abs(mask_right(1,kk) - mask_center(1,kk))
               rad_left_y  = abs(mask_left(2,kk)  - mask_center(2,kk))
               rad_right_y = abs(mask_right(2,kk) - mask_center(2,kk))
               rad_left_z  = abs(mask_left(3,kk)  - mask_center(3,kk))
               rad_right_z = abs(mask_right(3,kk) - mask_center(3,kk))
               !
               mask_radius(1,kk) = min(rad_left_x, rad_right_x)
               mask_radius(2,kk) = min(rad_left_y, rad_right_y)
               mask_radius(3,kk) = min(rad_left_z, rad_right_z)
               !
               if (ltrace3) then
                    do i = 1, 3
                        write(*,98)myid,' radius in each direction =',
     *              mask_radius(i,kk)
                     end do
               end if
               !mask_radius(1,kk) = 0.75d0*min(mask_radius(1,kk), 
               mask_radius(1,kk) = 0.5d0*min(mask_radius(1,kk), 
     &                            mask_radius(2,kk),mask_radius(3,kk))
               call MPI_BCAST(mask_radius(1,kk), 1, 
     *                  MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
               if(ltrace2)write(*,98)myid,'glb radius',mask_radius(1,kk)
               !
            else
               !
               if(ltrace2)write(*,99)myid,'No masking   for region: ',kk
               !
            end if
         end do  
         !
         ! Restore num_masks to proper value:
         num_masks = tmp_num_masks   
         do kk = 1, max_num_masks
            bh_true(kk) = tmp_bh_true(kk)
         end do
         call mijans_horizon(level)
      else if (abs(findhorizon).eq.6) then
         ! Do nothing -- see also findhorizon = 0
      else 
         !........................................................
         !      U N K N O W N    O P T I O N
         !........................................................
         if (myid.eq.master) write(*,*) 'Unknown value of findhorizon',
     *        findhorizon
         return
      end if

      !-------------------------------------------------------------
      ! Neutrino transport begins here
      !-------------------------------------------------------------
      if ( findhorizon .lt. 0 ) then
        call gpu_neutrino(level)
      end if
   
      !
      ! For debugging only:
      ! Copies hole #1 info into hole #2 but with the x-coordinate flipped:
      !
      if (forcesymmetric) then
           write(*,99) myid, 'level_bhmask_local: Forcing symmetric:'
           bh_true(2)       =  bh_true(1)
           mask_center(1,2) = -mask_center(1,1)
           mask_center(2,2) =  mask_center(2,1)
           mask_center(3,2) =  mask_center(3,1)
           mask_radius(1,2) =  mask_radius(1,1)
           !
           mask_bbox(1,2) = -mask_bbox(2,1)
           mask_bbox(2,2) = -mask_bbox(1,1)
           do i = 3,6
              mask_bbox(i,2) = mask_bbox(i,1)
           end do
      end if

      !
      ! Output the location of the masks (however located):
      !    (at this point in the code,
      !     all processors should have the same values)
      !
      if (ltrace2.or. (ltrace2.and. myid.eq.master) ) then
        do k = 1, max_num_masks
           if (bh_true(k)) then
              write(*,99) myid, 'level_bhmask_local: *****'	   
              write(*,99) myid, 'level_bhmask_local: Hole #: ',k
              do i = 1, 3
               write(*,98)myid,'    mask_center(i,k) =',mask_center(i,k)
              end do
              write(*,98)myid,'    mask_radius(1,k) =',mask_radius(1,k)
              write(*,98)myid,'      . . . . . '
              do i = 1, 6, 2
                 write(*,98)myid,'  min ---> max: ',
     *                                   mask_bbox(i,k),mask_bbox(i+1,k)
              end do
           end if
        end do
      end if


      !
      ! Finder finds a cube, but excision region is stored
      ! as a rectangle to allow for nice motion:
      !
      num_masks        = 0
      any_mask_changed = .false.
      do kk = 1, max_num_masks
         if (bh_true(kk)) then
            num_masks        = num_masks + 1
            mask_changed = .false.
            do i = 1, 6, 2
                oldmincoord = mask_coords(i,  kk)
                oldmaxcoord = mask_coords(i+1,kk)
                !
                ! Using new style:
                newmincoord = mask_bbox(i,  kk)
                newmaxcoord = mask_bbox(i+1,kk)
                !
                if ( noprevioushole(kk) ) then
                   ! This is a new region, no need to limit:
                   mask_coords(i,  kk) = newmincoord
                   mask_coords(i+1,kk) = newmaxcoord
                   mask_changed = .true.
                   if (ltrace.and.myid.eq.master)
     *                write(*,95)myid,'No prev. hole, so no limiting',kk
                else
                   !
                   if ( forcedmerge .and. LIMITSHRINKAGE
     *                  .and. kk.eq.1) then
                      !
                      ! Check if region is shrinking too much
                      ! without restricting its ability to move
                      !
                      newradius = newmaxcoord - newmincoord
                      if (.true. .and. myid.eq.master)
     *                write(*,95)myid,' Checking for shrinkage i/kk/r:',
     *                   i,kk,newradius
                      if (newradius .lt. MIN_AH_RADIUS) then
                         if (.true. .and. myid.eq.master)
     *                   write(*,95)myid,' Preventing shrinkagei/kk/r:',
     *                   i,kk,newradius
                         if (.true. .and. myid.eq.master)
     *                   write(*,98)myid,'  Before: newmin/maxcoord: ',
     *                   newmincoord,newmaxcoord
                         !
                         if (newmincoord .gt. oldmincoord .and.
     *                       newmaxcoord .lt. oldmaxcoord ) then
                            newmincoord = 0.5*(newmaxcoord+newmincoord)
     *                                    - MIN_AH_RADIUS
                            newmaxcoord = 0.5*(newmaxcoord+newmincoord)
     *                                    + MIN_AH_RADIUS
                            if (.true. .and. myid.eq.master)
     *                         write(*,95)myid,' Adjusting  both coords'
                         else if (newmincoord .gt. oldmincoord) then
                            newmincoord = newmaxcoord - 2.*MIN_AH_RADIUS
                            if (.true. .and. myid.eq.master)
     *                         write(*,95)myid,' Adjusting mincoord'
                         else 
                            newmaxcoord = newmincoord + 2.*MIN_AH_RADIUS
                            if (.true. .and. myid.eq.master)
     *                         write(*,95)myid,' Adjusting maxcoord'
                         end if
                         if (.true. .and. myid.eq.master)
     *                   write(*,98)myid,'  After: newmin/maxcoord: ',
     *                   newmincoord,newmaxcoord
                         !
                      end if
                   end if
                   !
                   ! Make sure mask does not change too quickly
                   if (ltrace.and.myid.eq.master)
     *                write(*,95)myid,'Checking if limits are needed',kk
                   if (newmincoord .gt. mask_coords(i,  kk)+hi) then
                      newmincoord =  mask_coords(i,kk) + hi
                      if (ltrace.and.myid.eq.master)
     *                write(*,95)myid,'ALimiting mask motion: ',i,kk,
     *                                newmincoord,mask_coords(i,kk),hi
                   else if(newmincoord.lt.mask_coords(i,  kk)-2*hi) then
                      newmincoord = mask_coords(i,kk) - 2.*hi
                      if (ltrace.and.myid.eq.master)
     *                write(*,95)myid,'ALimiting mask motion: ',i,kk,
     *                                newmincoord,mask_coords(i,kk),hi
                   else
                      !mask_coords(i,  kk) = newmincoord
                      if (ltrace.and.myid.eq.master)
     *                write(*,95)myid,'ANot limitmask motion: ',i,kk,
     *                                newmincoord,mask_coords(i,kk),hi
                   end if
                   if (newmaxcoord .lt. mask_coords(i+1,kk)-hi) then
                      newmaxcoord = mask_coords(i+1,kk) - hi
                      if (ltrace.and.myid.eq.master)
     *                write(*,95)myid,'BLimiting mask motion: ',i,kk,
     *                                newmaxcoord,mask_coords(i+1,kk),hi
                   else if(newmaxcoord.gt.mask_coords(i+1,kk)+2*hi) then
                      newmaxcoord =  mask_coords(i+1,kk) + 2.*hi
                      if (ltrace.and.myid.eq.master)
     *                write(*,95)myid,'BLimiting mask motion: ',i,kk,
     *                                newmaxcoord,mask_coords(i+1,kk),hi
                   else
                      !mask_coords(i+1,kk) = newmaxcoord
                      if (ltrace.and.myid.eq.master)
     *                write(*,95)myid,'BNot limitmask motion: ',i,kk,
     *                                newmaxcoord,mask_coords(i+1,kk),hi
                   end if
                   !
                   ! Set new max region now that limits have been placed:
                   !
                   mask_coords(i,  kk) = newmincoord
                   mask_coords(i+1,kk) = newmaxcoord
                   !
                   mask_changed = mask_changed .or.
     *                              (oldmincoord.ne.mask_coords(i,  kk))
     *                          .or.(oldmaxcoord.ne.mask_coords(i+1,kk))
                end if
            end do
            if (mask_changed .and. (myid .eq. master.or.ltrace2)) then
              write(*,99) myid,'level_bhmask_local: Changed Hole #: ',kk
              do i = 1, 6, 2
                 write(*,98)myid,'       mask_coords(i,kk)=',
     *                             mask_coords(i,kk),mask_coords(i+1,kk)
              end do
            end if
            if (mask_changed.and.myid.eq.master) then
                open(unit=12,file='maskcoords.dat',STATUS='UNKNOWN')
                write(12,"(I4,4F11.5)")  kk,time,
     *                       0.5*(mask_coords(1,kk)+mask_coords(2,kk)),
     *                       0.5*(mask_coords(3,kk)+mask_coords(4,kk)),
     *                       0.5*(mask_coords(5,kk)+mask_coords(6,kk))
                open(unit=13,file='maskbbox.dat',STATUS='UNKNOWN')
!               write(13,"(I4,4F9.5)")  kk,time,mask_center(1,kk),
!    *             mask_center(2,kk),mask_center(3,kk),mask_radius(1,kk)
                write(13,"(I4,7F10.5)")  kk,time,
     *                              mask_bbox(1,kk),mask_bbox(2,kk),
     *                              mask_bbox(3,kk),mask_bbox(4,kk),
     *                              mask_bbox(5,kk),mask_bbox(6,kk)
            end if
            any_mask_changed = any_mask_changed .or. mask_changed
         end if
      end do


      ! ........................................................
      ! Try to force the holes to MERGE if close enough?
      ! ........................................................
      if (forceclosemerge .and. bh_true(1) .and. bh_true(2)) then
           overlap    = 1.d0
           BHdistance = 0.d0
           do i = 1,6,2
              ! distance between centers:
              BHdistance = BHdistance +
     *                 0.25d0*( mask_coords(i+1,2)+mask_coords(i,2)
     *                         -mask_coords(i+1,1)-mask_coords(i,1))**2
              ! Normalized fractional overlap in each direction:
	      oversing = min( mask_coords(i+1,2)-mask_coords(i,1),
     *                       mask_coords(i+1,1)-mask_coords(i,2)  )

c              overlap  = overlap *
c     *                  min( mask_coords(i+1,2)-mask_coords(i,1),
c     *                       mask_coords(i+1,1)-mask_coords(i,2)  )
c     *                  /  
c     *                  max( mask_coords(i+1,1)-mask_coords(i,1),
c     *                       mask_coords(i+1,2)-mask_coords(i,2)  )
               !
               ! If overlap is negative, it means there is no overlap:
               ! unless it is tiny in which case we try to go on
       if (oversing.lt.0.and.abs(overlap).lt.small_delta*5.) oversing=0.01d0


	       overlap = overlap * oversing


               if (overlap .lt. 0) overlap = 0.d0
               if (ltrace2) write(*,96)myid, ' overlap:   ',i,overlap
               if (ltrace2) write(*,96)myid, ' BHdistance:',i,BHdistance
           end do
           BHdistance = sqrt(BHdistance)
           if (ltrace2.and. myid.eq.master) then
           write(*,98)myid,'level_bhmask_local:overlap:   ',overlap
           write(*,98)myid,'level_bhmask_local:BHdistance:',BHdistance
           end if
           if (     overlap .gt. MERGETHRESH
     *                      .or.
     *              BHdistance .lt. MERGETHRESHDistance) then
              !
              if(ltrace)write(*,99)myid,'level_bhmask_local:Force merge'
              !
              num_masks        =  1
              mask_changed     = .true.
              any_mask_changed = .true.
              !
              forcedmerge      = .true.
              transition_time  = time
              !
              ! Get rid of second masked region:
              bh_true(2) = .false.
              !
              ! Set first masked region to Union of both:
              !
              mask_radius(1,1) = LARGENUMBER
              do i = 1, 6, 2
                 mask_coords(i,1)   = min( mask_coords(i,  1),  
     *                                     mask_coords(i,  2)  )
                 mask_coords(i+1,1) = max( mask_coords(i+1,1),
     *                                     mask_coords(i+1,2)  )
                 !
                 mask_bbox(i,  1)   = mask_coords(i,  1)
                 mask_bbox(i+1,1)   = mask_coords(i+1,1)
                 !
                 mask_center(i,1)=.5d0*(mask_bbox(i+1,1)+mask_bbox(i,1))
                 mask_radius(1,1)= min( mask_radius(1,1), 
     *                         0.5d0*(mask_bbox(i+1,1)-mask_bbox(i,1)) )
                 if (ltrace2) then
                    write(*,96)myid,'New bounds: ',i,
     *                                   mask_bbox(i,1),mask_bbox(i+1,1)
                    !write(*,98)myid,'Ctr: ',mask_center(i,1)
                 end if
              end do 
              !
              if (myid.eq.master) then
                 write(*,99)myid,'level_bhmask_local: Merged Hole #: ',1
                 do i = 1, 6, 2
                    write(*,98)myid,'       mask_coords(i,1)=',
     *                             mask_coords(i,1),mask_coords(i+1,1)
                 end do
              end if
              !
           end if
      end if

      !
      ! Only need to mask if:
      !     (1) the mask has changed
      !  or
      !     (2) the level changed, in which case we have to mask on this level
      !
      if(ltrace.and.myid.eq.master)
     *      write(*,99)myid,'level_bhmask_local:num_masks',num_masks
      if ( num_masks .gt. 0    .and.
     *     (any_mask_changed .or. levelchanged) ) then
         if (ltrace.and.myid.eq.master)
     *   write(*,99)myid,'level_bhmask_local: calling grid_bhmask()'
         !
         ! Implement the mask:
         !
         gi = level_return_start(level)
  10     if (grid_return_existence(gi)) then
            if (grid_is_local(gi)) then
               if (ltrace2) write(*,99)myid,'level_bhmask_local: gi=',gi
               call grid_bhmask(gi)
            end if
            gi = grid_return_sibling(gi)
            goto 10
         end if
      end if

      !
      ! If we had masked regions before, and now we are not
      ! finding them, then let us return (but need to retain
      ! the masks, they cannot just disappear):
      !     (except that two masks can merge)
      !
      num_masks = 0
      do kk = 1, max_num_masks
         if (bh_true(kk)) then
            num_masks = num_masks + 1
         else if ( .not.    noprevioushole(kk)  .and.
     *             .not.(kk.eq.2.and.forcedmerge)      ) then
         !else if (.not.noprevioushole(kk)) then
            if (myid.eq.master)
     *         write(*,99)myid,'level_bhmask_local: Hole vanished!',kk
            bh_true(kk) = .true.
            num_masks = num_masks + 1
         end if
      end do

      gridnum=19
      if(ltrace.and.grid_is_local(gridnum)) then
            call field_out3d(q(gfunc_pointer(outfnum,gridnum)),
     *             gr_t(gridnum),'BHpos',
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             gr_nx(gridnum),gr_ny(gridnum),gr_nz(gridnum),myid)
            call field_out3d(q(gfunc_pointer(outfnum+1,gridnum)),
     *             gr_t(gridnum),'BHposADV',
     *             gr_minx(gridnum),gr_maxx(gridnum),
     *             gr_miny(gridnum),gr_maxy(gridnum),
     *             gr_minz(gridnum),gr_maxz(gridnum),
     *             gr_nx(gridnum),gr_ny(gridnum),gr_nz(gridnum),myid)
      end if
      if (ltrace2) write(*,99) myid, 'level_bhmask_local: Finished'

  95  format('[',I3,'] ',A,2I5,3F18.10)
  96  format('[',I3,'] ',A,I5,3F18.10)
  98  format('[',I3,'] ',A,3F18.10)
  99  format('[',I3,'] ',A,3I5)
      return 
      end    ! END: level_bhmask_local

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_init_local:                                                         cc
cc                       Each processor initializes each grid local to        cc
cc                  it.                                                       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_init_local(level)
      implicit     none
      integer      level
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'glob.inc'
      include     'grid_methods.inc'
      integer      gi

      logical      ltrace
      parameter (  ltrace = .false. )

      !
      gi = level_return_start(level)
      !
  10  if (grid_return_existence(gi)) then
         if (grid_is_local(gi)) then
            if (ltrace) write(*,99) myid, 'level_init_local: gi=',gi
            call grid_initialize(gi)
         end if
         gi = grid_return_sibling(gi)
         goto 10
      end if

      if (ltrace) call level_test_local()
      if (ltrace) write(*,99) myid, 'level_init_local: Finished'

  99  format('[',I3,'] ',A,3I5)
  98  format('[',I3,'] ',A,3F18.10)

      return 
      end    ! END: level_init_local


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_resettime_local:                                                    cc
cc                  Set t=0 for each grid local to proc in the level.         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_resettime_local(level)
      implicit     none
      integer      level
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      include     'glob.inc'    ! for lev_count
      integer      gi

      logical      ltrace
      parameter (  ltrace = .false. )

      !
      gi = level_return_start(level)
      !
  10  if (grid_return_existence(gi)) then
         !if (grid_is_local(gi)) then
            if (ltrace)write(*,99) myid, 'level_resettime_local: gi=',gi
            call grid_reset_time(gi)
         !end if
         gi = grid_return_sibling(gi)
         goto 10
      end if

      if(ltrace)write(*,99)myid,'level_resettime_local:Resetg lev_count'
      if(ltrace)write(*,99) myid, 'level_resettime_local: old:',
     *        lev_count(level)
      if (level.ge.0 .and. level .lt. maxlev) then
         lev_count(level) = 0
      else
         write(*,99)myid,'level_resettime: Something wrong:',level
      end if

      if (ltrace) write(*,99) myid, 'level_resettime_local: Finished'

  99  format('[',I3,'] ',A,3I5)
  98  format('[',I3,'] ',A,3F18.10)

      end    ! END: level_resettime_local

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_cgctransfer_local:                                                  cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_cgctransfer_local(level)
      implicit     none
      integer      level
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      gi, gj, gi_owner, gj_owner
      integer      round, numrounds, proc
      integer      round_robin
      external     round_robin
      logical      gi_local, gj_local
      logical      ltrace
      parameter (  ltrace = .false. )

      if (ltrace) then
         write(*,99) myid, 'level_cgctransfer_local: level = ',level
      end if

      !
      ! We are processor "myid" and we need to check with all other
      ! processors to see if we need to sync with them. We want to do
      ! so in a way that avoids having a processor wait for the other
      ! processor to finish with someone else. So we use a round robin.
      ! 
      ! For each "round" of the round robin, we determine which processor
      ! with which to sync.

      if (mod(numprocs,2).eq.0) then
         ! EVEN number of processors:
         numrounds = numprocs - 1
      else
         ! ODD  number of processors:
         numrounds = numprocs
      end if
      ! Add a round to sync pairs of grids both owned by this process:
      numrounds = numrounds + 1

      !
      ! Loop over the rounds of the round robin:
      !
      do round = 1, numrounds
        if(ltrace)write(*,99)myid,'level_cgctransfer_local:round=',round
        !
        ! Determine which processor with which to sync:
        !
        if (round .ne. numrounds) then
           proc = round_robin(round, myid, numprocs)
        else
           proc = myid
        end if
        !
        ! Make sure we don't have a "bye" (only happens for odd numprocs):
        !
        if (proc .ge. 0) then
           if (ltrace) then
              write(*,99)myid,'level_cgctransfer_local: w/proc=',proc
           end if
           !
           ! Loop for all grids (gi) on level. For each gi,
           ! consider pairings of gi and gj where occurs in the
           ! hierarchy after gi (this avoids duplicate pairings),
           ! and sync those pairings (gi,gj):
           !
           gi = level_return_start(level)
           !
 10        if (grid_return_existence(gi)) then
              if (ltrace) write(*,99) myid, '      gi = ',gi
              gi_owner = grid_return_owner(gi)
              gi_local = grid_is_local(gi)
              ! If neither us nor processor "proc" own gi, then skip:
              if (.not. gi_local .and. gi_owner .ne. proc) goto 40
              !
              !gj = grid_return_sibling(gi)
              ! Loop over parent level:
              gj = level_return_start(level-1)
 20           if (grid_return_existence(gj)) then
                 !
                 ! Should we skip this pairing?
                 gj_owner = grid_return_owner(gj)
                 gj_local = grid_is_local(gj)
                 if (.not. gj_local .and. gj_owner .ne. proc) goto 30
                 if (gi_local.and.gj_local .and. (myid.ne.proc)) goto 30
                 ! Does pair intersect?
                 if (ltrace) write(*,99) myid, '      w/ gj = ',gj
                 if ( grid_intersect(gi,gj) ) then
                    if (gj_local) then
                       if (gi_local) then
                          if (ltrace) write(*,99) myid, '   gi/gj local'
                          !call grid_inject(gi, gj)
                          call grid_cgctransfer(gi, gj)
                       else
                          if (ltrace) write(*,99) myid, '   gj   local'
                          !call grid_inject_master(gi, gj)
                          call grid_cgctransfer_master(gi, gj)
                       end if
                    else
                       if (gi_local) then
                          if (ltrace) write(*,99) myid, '   gi   local'
                          !call grid_inject_slave(gi, gj)
                          call grid_cgctransfer_slave(gi, gj)
                       else
                          if (ltrace) write(*,99) myid, '   none local'
                       end if
                    end if
                 else
                          if (ltrace) write(*,99) myid, '   No intersec'
                 end if
                 ! Look to next sibling
  30             gj = grid_return_sibling(gj)
                 goto 20
              end if
              !
              ! No more siblings w/ which to sync gi,
              ! so, repeat process for next grid on level:
              !
  40          gi = grid_return_sibling(gi)
              goto 10
           end if
           !
        end if
      end do

      if (ltrace) write(*,99) myid, 'level_cgctransfer_local: Finished.'

  99  format('[',I3,'] ',A,3I5)

      return 
      end    ! END: level_cgctransfer_local

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_psttransfer_local:                                                  cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_psttransfer_local(level)
      implicit     none
      integer      level
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      gi, gj, gi_owner, gj_owner
      integer      round, numrounds, proc
      integer      round_robin
      external     round_robin
      logical      gi_local, gj_local
      logical      ltrace
      parameter (  ltrace = .false. )

      if (ltrace) then
         write(*,99) myid, 'level_psttransfer_local: level = ',level
      end if

      !
      ! We are processor "myid" and we need to check with all other
      ! processors to see if we need to sync with them. We want to do
      ! so in a way that avoids having a processor wait for the other
      ! processor to finish with someone else. So we use a round robin.
      ! 
      ! For each "round" of the round robin, we determine which processor
      ! with which to sync.

      if (mod(numprocs,2).eq.0) then
         ! EVEN number of processors:
         numrounds = numprocs - 1
      else
         ! ODD  number of processors:
         numrounds = numprocs
      end if
      ! Add a round to sync pairs of grids both owned by this process:
      numrounds = numrounds + 1

      !
      ! Loop over the rounds of the round robin:
      !
      do round = 1, numrounds
        if(ltrace)write(*,99)myid,'level_psttransfer_local:round=',round
        !
        ! Determine which processor with which to sync:
        !
        if (round .ne. numrounds) then
           proc = round_robin(round, myid, numprocs)
        else
           proc = myid
        end if
        !
        ! Make sure we don't have a "bye" (only happens for odd numprocs):
        !
        if (proc .ge. 0) then
           if (ltrace) then
              write(*,99)myid,'level_psttransfer_local: w/proc=',proc
           end if
           !
           ! Loop for all grids (gi) on level. For each gi,
           ! consider pairings of gi and gj where occurs in the
           ! hierarchy after gi (this avoids duplicate pairings),
           ! and sync those pairings (gi,gj):
           !
           gi = level_return_start(level)
           !
 10        if (grid_return_existence(gi)) then
              if (ltrace) write(*,99) myid, '      gi = ',gi
              gi_owner = grid_return_owner(gi)
              gi_local = grid_is_local(gi)
              ! If neither us nor processor "proc" own gi, then skip:
              if (.not. gi_local .and. gi_owner .ne. proc) goto 40
              !
              !gj = grid_return_sibling(gi)
              ! Loop over parent level:
              gj = level_return_start(level-1)
 20           if (grid_return_existence(gj)) then
                 !
                 ! Should we skip this pairing?
                 gj_owner = grid_return_owner(gj)
                 gj_local = grid_is_local(gj)
                 if (.not. gj_local .and. gj_owner .ne. proc) goto 30
                 if (gi_local.and.gj_local .and. (myid.ne.proc)) goto 30
                 ! Does pair intersect?
                 if (ltrace) write(*,99) myid, '      w/ gj = ',gj
                 if ( grid_intersect(gi,gj) ) then
                    if (gj_local) then
                       if (gi_local) then
                          if (ltrace) write(*,99) myid, '   gi/gj local'
                          !call grid_inject(gi, gj)
                          call grid_psttransfer(gi, gj)
                       else
                          if (ltrace) write(*,99) myid, '   gj   local'
                          !call grid_inject_master(gi, gj)
                          call grid_psttransfer_master(gi, gj)
                       end if
                    else
                       if (gi_local) then
                          if (ltrace) write(*,99) myid, '   gi   local'
                          !call grid_inject_slave(gi, gj)
                          call grid_psttransfer_slave(gi, gj)
                       else
                          if (ltrace) write(*,99) myid, '   none local'
                       end if
                    end if
                 else
                          if (ltrace) write(*,99) myid, '   No intersec'
                 end if
                 ! Look to next sibling
  30             gj = grid_return_sibling(gj)
                 goto 20
              end if
              !
              ! No more siblings w/ which to sync gi,
              ! so, repeat process for next grid on level:
              !
  40          gi = grid_return_sibling(gi)
              goto 10
           end if

        end if
      end do

      if (ltrace) write(*,99) myid, 'level_psttransfer_local: Finished.'

  99  format('[',I3,'] ',A,3I5)
      return 
      end    ! END: level_psttransfer_local

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_finest_containing:                                                  cc
cc              Finds the level and grid at the finest level which            cc
cc              contains the passed in point (x,y,z).                         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_finest_containing(x,y,z,lev,grid,owner)
      implicit    none
      real(kind=8)x,y,z
      integer     lev, grid, owner
      include    'glob.inc'
      include    'grid_methods.inc'
      logical     pointfound

      !
      ! Start on finest level that exists now;
      !
      lev = level_return_finest(Levelp, maxlev)

  9   grid = level_return_start(lev)
 10   if (grid_return_existence(grid)) then
         !
         ! Is point on this grid?
         !
         pointfound = grid_contains_pt(grid,x,y,z)
         !
         if (.not. pointfound) then
            grid = grid_return_sibling(grid)
            goto 10
         end if
      end if

      if (.not.pointfound)  then
         lev = lev - 1
         if (lev.ge.0) goto 9
      end if 


      if (.not.pointfound)  then
         lev   = -1
         grid  = -1
         owner = -1
      else
         owner = grid_return_owner(grid)
      end if

      return
      end       ! END: level_finest_containing

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  grid_minfind:                                                             cc
cc                 Find the minimum of a given field on this grid.            cc
cc                 If less than the input minimum, update the coordinates.    cc
cc                 Variable "num" determines whether to search x<0 or x>0.    cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_minfind(grid,fieldnum,mymin,coords,num)
      implicit     none
      integer      grid, fieldnum, num
      real(kind=8) mymin, coords(3)
      include     'grid.inc'
      integer      nx, ny, nz, i, j, k, l, index
      real(kind=8) tmpmin, tmpcoords(3), x, value
      include     'largesmall.inc'
      !real(kind=8)   LARGENUMBER
      !parameter    ( LARGENUMBER = 9.d98 )

      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )

      !call load_pointers(grid)

      nx = gr_nx(grid)
      ny = gr_ny(grid)
      nz = gr_nz(grid)

      tmpmin = LARGENUMBER

      if (ltrace) then
         write(*,*) 'grid_minfind:             grid: ',grid
         write(*,*) 'grid_minfind:            mymin: ',mymin
         write(*,*) 'grid_minfind: Outputting field: ',fieldnum
         write(*,*) 'grid_minfind: nx/y/z:           ',nx,ny,nz
         write(*,*) 'grid_minfind: min/max: ',
     *             gr_minx(grid),gr_maxx(grid),
     *             gr_miny(grid),gr_maxy(grid),
     *             gr_minz(grid),gr_maxz(grid)
         call field_out3d(q(gfunc_pointer(fieldnum,grid)),gr_t(grid),
     *             'grid_minfind_field',
     *             gr_minx(grid),gr_maxx(grid),
     *             gr_miny(grid),gr_maxy(grid),
     *             gr_minz(grid),gr_maxz(grid),
     *             nx,ny,nz,0)
      end if

      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         index = (k-1)*nx*ny + (j-1)*nx + (i-1)
         !
         value = q(gfunc_pointer(fieldnum,grid)+index)
         x     = q(gr_x(grid)+i-1)
         !
         if (    ( (num.eq.1) .and. (x.le.0) )
     *                  .or.
     *           ( (num.eq.2) .and. (x.ge.0) )    ) then
            if (value .lt. tmpmin) then
               tmpmin       = value
               tmpcoords(1) = x
               tmpcoords(2) = q(gr_y(grid)+j-1)
               tmpcoords(3) = q(gr_z(grid)+k-1)
               if (ltrace2) then
                  write(*,*) 'grid_minfind:     New min found:',tmpmin
                  write(*,*) '               tmpcoords(1):',tmpcoords(1)
                  write(*,*) '               tmpcoords(2):',tmpcoords(2)
                  write(*,*) '               tmpcoords(3):',tmpcoords(3)
               end if
            end if
         end if
      end do
      end do
      end do

      if (tmpmin .lt. mymin) then
         mymin = tmpmin
         do l = 1, 3
            coords(l) = tmpcoords(l)
         end do
         if (ltrace) then
            write(*,*) 'grid_minfind:     New min found:',mymin
            write(*,*) '               coords(1):',coords(1)
            write(*,*) '               coords(2):',coords(2)
            write(*,*) '               coords(3):',coords(3)
         end if
      else
         if (ltrace) write(*,*) 'grid_minfind:     No new min found'
      end if

      return
      end    ! END: grid_minfind

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_analysis_local:                                                     cc
cc                        Carryout various analysis functions on              cc
cc                        the fields....not necessary for evolution.          cc
cc                                                                            cc
cc        == 0  [default] no analysis                                         cc
cc        == 1  spatial integral                                              cc
cc        == 2  l2norm                                                        cc
cc        == 4  global min                                                    cc
cc        == 8  global max                                                    cc
cc        ==16  global min of abs()                                           cc
cc        ==32  global max of abs()                                           cc
cc        ==64  first and second moments (com and quadrupole)                 cc
cc        ==128 central value
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_analysis_local(level)
      implicit     none
      integer      level
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      i, option, numfields
      integer      MIN,     MAX,     ABSMIN,     ABSMAX
      parameter (  MIN = 0, MAX = 1, ABSMIN = 2, ABSMAX = 3 )
      integer      INTEGRATE,     L2NORM
      parameter (  INTEGRATE = 2, L2NORM = 3 )
      logical      masksetyet

      logical      ltrace
      parameter (  ltrace = .false. )

      if (ltrace)write(*,99) myid, 'level_analysis_local: level= ',level

      masksetyet = .false.
      numfields  = grid_return_numgfuncs()
      !
      do i = 1, numfields
         option = grid_return_ana(i)
         if (option .ne. 0) then
            !
            if (.not.masksetyet) then
               ! Set mask if it has not been set yet 
               !    (if no fields require analysis, no need to set):
               call level_set_anamask(level)
               masksetyet = .true.
            end if
            !
            if(ltrace) then
                write(*,99)myid,'level_analysis_local:Analyze field i',i
                write(*,99)myid,'level_analysis_local: Option: ',option
            end if
            !
            if (mod(option,2) .eq. 1) then
               ! If odd, take the integral:
               if(ltrace)write(*,99)myid,'         ...Integrate'
               call level_comp_integral(level,i,INTEGRATE)
            end if
            if ((option.ge.    2 .and. option.lt. 2* 2) .or.
     *          (option.ge. 3* 2 .and. option.lt. 4* 2) .or.
     *          (option.ge. 5* 2 .and. option.lt. 6* 2) .or.
     *          (option.ge. 7* 2 .and. option.lt. 8* 2) .or.
     *          (option.ge. 9* 2 .and. option.lt.10* 2) .or.
     *          (option.ge.11* 2 .and. option.lt.12* 2) .or.
     *          (option.ge.13* 2 .and. option.lt.14* 2) .or. 
     *          (option.ge.15* 2 .and. option.lt.16* 2) .or. 
     *          (option.ge.17* 2 .and. option.lt.18* 2) .or. 
     *          (option.ge.19* 2 .and. option.lt.20* 2) .or. 
     *          (option.ge.21* 2 .and. option.lt.22* 2) .or. 
     *          (option.ge.23* 2 .and. option.lt.24* 2) .or. 
     *          (option.ge.25* 2 .and. option.lt.26* 2) .or. 
     *          (option.ge.27* 2 .and. option.lt.28* 2) .or. 
     *          (option.ge.29* 2 .and. option.lt.30* 2) .or. 
     *          (option.ge.31* 2 .and. option.lt.32* 2) .or.
     *          (option.ge.32* 2 .and. option.lt.34* 2) ) then
               if(ltrace)write(*,99)myid,'         ...l2norm'
               call level_comp_integral(level,i,L2NORM)
            end if
            if ((option.ge.    4 .and. option.lt. 2* 4) .or.
     *          (option.ge. 3* 4 .and. option.lt. 4* 4) .or.
     *          (option.ge. 5* 4 .and. option.lt. 6* 4) .or.
     *          (option.ge. 7* 4 .and. option.lt. 8* 4) .or.
     *          (option.ge. 9* 4 .and. option.lt.10* 4) .or.
     *          (option.ge.11* 4 .and. option.lt.12* 4) .or.
     *          (option.ge.13* 4 .and. option.lt.15* 4) .or. 
     *          (option.ge.15* 4 .and. option.lt.16* 4) .or.
     *          (option.ge.17* 4 .and. option.lt.18* 4) ) then
               if(ltrace)write(*,99)myid,'         ...global min'
               call level_comp_extremum(level,i,MIN)
            end if
            if ((option.ge.   8 .and. option.lt. 2* 8) .or.
     *          (option.ge.3* 8 .and. option.lt. 4* 8) .or.
     *          (option.ge.5* 8 .and. option.lt. 6* 8) .or.
     *          (option.ge.7* 8 .and. option.lt. 8* 8) .or.
     *          (option.ge.9* 8 .and. option.lt.10* 8) ) then
               if(ltrace)write(*,99)myid,'         ...global max'
               call level_comp_extremum(level,i,MAX)
            end if
            if ((option.ge.  16 .and. option.lt.2*16) .or.
     *          (option.ge.3*16 .and. option.lt.4*16) .or.
     *          (option.ge.5*16 .and. option.lt.6*16) ) then
               if(ltrace)write(*,99)myid,'         ...global absmin'
               call level_comp_extremum(level,i,ABSMIN)
            end if
            if ((option.ge.  32 .and. option.lt.2*32) .or.
     *          (option.ge.3*32 .and. option.lt.4*32) ) then
               if(ltrace)write(*,99)myid,'         ...global absmax'
               call level_comp_extremum(level,i,ABSMAX)
            end if
            if ((option.ge.  64 .and. option.lt.2*64) ) then
               if(ltrace)write(*,99)myid,'         ...moments'
               call level_comp_moments(level,i)
            end if
            if ((option.ge. 128 .and. option.lt.2*128) ) then
               if(ltrace)write(*,99)myid,'         ...central'
               call level_central_value(level,i)
            end if
            !
         end if
      end do

      if (ltrace) write(*,99) myid, 'level_analysis_local: Finished'

  98  format('[',I3,'] ',A,F15.8)
  99  format('[',I3,'] ',A,3I5)
      return 
      end    ! END: level_analysis_local

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_set_anamask:                                                        cc
cc                      Set a temporary mask for use in computing             cc
cc                      the analysis_XXX functions. Not related to            cc
cc                      either the chr() array or the black hole mask.        cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_set_anamask(level)
      implicit     none
      integer      level
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      lev, gi, i

      logical      ltrace
      parameter (  ltrace  = .false. )
      logical      ltrace2
      parameter (  ltrace2 = .false. )

      if (ltrace)write(*,99) myid, 'level_set_anamask:   level: ',level

      !
      lev     = level
  10  if (level_return_existence(lev)) then
         gi = level_return_start(lev)
  15     if (grid_return_existence(gi)) then
            if (grid_is_local(gi)) then
               if(ltrace2)write(*,99)myid,'call grid_set_anamask: ',gi
               call grid_set_anamask(gi)
            end if
            gi = grid_return_sibling(gi)
            goto 15
         end if
         lev = lev + 1
         goto 10
      end if
      !
      if (ltrace) then
         write(*,99)myid,'level_set_anamask: Done.'
      end if

  97  format('[',I3,'] ',A,I4,F15.8)
  98  format('[',I3,'] ',A,F15.8)
  99  format('[',I3,'] ',A,3I5)

      return 
      end    ! END: level_set_anamask

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_comp_extremum:                                                      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_comp_extremum(level,field,minormax)
      implicit     none
      integer      level, field, minormax
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      MIN,     MAX,     ABSMIN,     ABSMAX
      parameter (  MIN = 0, MAX = 1, ABSMIN = 2, ABSMAX = 3 )
      integer      lev, gi, i
      integer      mpi_opt
      ! Processor ID of who owns the min or max:
      integer      procext
      real(kind=8) locext,  loccoords(3)
      real(kind=8) globext, globcoords(3)
      real(kind=8)   LARGENUMBER
      parameter    ( LARGENUMBER = 9.d98 )
      logical      double_equal
      external     double_equal
      ! Output variables:
      integer      mylen, gft_rc, rank, shape(3)
      parameter  ( rank = 3)
      data         shape / 1, 1, 1 /
      real(kind=8) time, bbox(6)
      character(128) fname,  procname, fieldnum
      ! Count the number of points over which l2norm is taken
      integer      tmppts, numpts, globnum
      integer      gft_out_bbox, mystringlength
      external     gft_out_bbox, mystringlength

      logical      ltrace
      parameter (  ltrace  = .false. )
      logical      ltrace2
      parameter (  ltrace2 = .false. )

      time = grid_return_time(level_return_start(level))
      if (ltrace) then
         write(*,99) myid, '----------------------------------'
         write(*,99) myid, 'level_comp_extremum: level: ',level
      end if

      !
      ! Find extremum in some field (say lapse):
      !    (1) Find local min/max on all levels as fine or finer than level
      !    (2) Find global min across all procs
      !
      if (minormax .eq. MAX .or. minormax .eq. ABSMAX) then
         locext  = -LARGENUMBER
         mpi_opt = MPI_MAX
      else
         locext  =  LARGENUMBER
         mpi_opt = MPI_MIN
      end if
      !
      lev     = level
  10  if (level_return_existence(lev)) then
         gi = level_return_start(lev)
  15     if (grid_return_existence(gi)) then
            if (grid_is_local(gi)) then
               if(ltrace2)write(*,99)myid,'call grid_extfind: ',gi
               call grid_extfind(gi, field, locext,loccoords, minormax)
            end if
            gi = grid_return_sibling(gi)
            goto 15
         end if
         lev = lev + 1
         goto 10
      end if
      !
      if (ltrace) then
         write(*,99)myid,'level_comp_extremum: Local min/max found'
         write(*,98)myid,'level_comp_extremum: locext:',locext
         do i = 1, 3
            write(*,98)myid,'level_comp_extremum: coords: ',loccoords(i)
         end do
      end if
      !
      ! Get the global answer:
      call MPI_Reduce(locext,   globext, 1, MPI_DOUBLE_PRECISION,
     *                mpi_opt,  master,  MPI_COMM_WORLD, ierr)
      if(ltrace.and.myid.eq.master)write(*,98)myid,' global: ',globext

      ! Tell all procs the answer:
      call MPI_BCAST(globext, 1, MPI_DOUBLE_PRECISION,
     &               master, MPI_COMM_WORLD, ierr)
      if(ltrace)write(*,98)myid,' globext: ',globext

      ! All procs have the global extremum, now figure out 
      ! if we have it:
      if (double_equal(locext,globext)) then
         ! We have it so put our processor id to all:
         call MPI_Reduce(myid,procext,1,MPI_INTEGER,
     *                         MPI_MIN,master,MPI_COMM_WORLD,ierr)
      else
         ! We do not have it, so just put out a big number:
         call MPI_Reduce(numprocs,procext,1,MPI_INTEGER,
     *                         MPI_MIN,master,MPI_COMM_WORLD,ierr)
      end if
      ! Let everyone know that procext will broadcast coords:
      call MPI_BCAST(procext, 1, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
      if (procext .lt. numprocs .and. procext.ge.0) then
         if(ltrace)write(*,99)myid,' procext has it: ',procext
         ! Everyone now gets the coords of the min from procext:
         call MPI_BCAST(loccoords,      3, MPI_DOUBLE_PRECISION,
     &                            procext, MPI_COMM_WORLD, ierr)
      else
         write(*,99)myid,' procext does not have a valid value'
         write(*,99)myid,' Something wrong'
         write(*,99)myid,' procext:  ',procext
         write(*,99)myid,' numprocs: ',numprocs
         write(*,98)myid,' locext:   ',locext
         write(*,98)myid,' globext:  ',globext
      end if

      do i = 1, 3
         bbox(2*i-1) = loccoords(i)
         bbox(2*i  ) = loccoords(i)
      end do
      !
      ! Output the min or max:
      !
      if (myid.eq.master) then
         call int2str(field,fieldnum)
         call grid_get_fname(field,fname)
         mylen           = mystringlength(fname)
         if (minormax .eq. MIN) then
            fname(mylen+1:) = '_min'
         else if (minormax .eq. ABSMIN) then
            fname(mylen+1:) = '_absmin'
         else if (minormax .eq. ABSMAX) then
            fname(mylen+1:) = '_absmax'
         else
            fname(mylen+1:) = '_max'
         end if
         if (ltrace) then
            write(*,99)myid,' procext:  ',procext
            write(*,99)myid,' numprocs: ',numprocs
            write(*,98)myid,' locext:   ',locext
            write(*,98)myid,' globext:  ',globext
            write(*,97)myid,' Coords:   ',1,bbox(1)
            write(*,97)myid,' Coords:   ',2,bbox(3)
            write(*,97)myid,' Coords:   ',3,bbox(5)
            write(*,*) myid,' fname: ',fname(1:13)
         end if
#ifdef AIX
         gft_rc = gft_out_bbox( trim(fname)/ /CHAR(0), %VAL(time),
     *                    shape, %VAL(rank), bbox, globext )
#else
         gft_rc = gft_out_bbox( fname, time, shape, rank, bbox, globext)
#endif
      end if

      if (ltrace) then
         write(*,98) myid, 'level_comp_extremum: Finished',globext
         do i = 1, 3
            write(*,97) myid, 'level_comp_extremum: ',i,loccoords(i)
         end do
         write(*,99) myid, '----------------------------------'
      end if

  97  format('[',I3,'] ',A,I4,F15.8)
  98  format('[',I3,'] ',A,F15.8)
  99  format('[',I3,'] ',A,3I5)
      return 
      end    ! END: level_comp_extremum

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_comp_integral:                                                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_comp_integral(level,field,intorl2norm)
      implicit     none
      integer      level, field, intorl2norm
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      INTEGRATE,     L2NORM
      parameter (  INTEGRATE = 2, L2NORM = 3 )
      integer      lev, gi, i, mylen, gft_rc
      real(kind=8) locint, globint, time, bbox(2)
      data         bbox / 0.d0, 0.d0 /
      character(128) fname,  procname, fieldnum
      ! Count the number of points over which l2norm is taken
      integer      tmppts, numpts, globnum
      integer      gft_out_bbox, mystringlength
      external     gft_out_bbox, mystringlength

      logical      ltrace
      parameter (  ltrace = .false. )

      if (ltrace)write(*,99) myid, 'level_comp_integral: level: ',level

      time = grid_return_time(level_return_start(level))
      !
      !  We need to:
      !     1) For each local grid on this level:
      !        a) mask off location of other grids on level
      !        b) mask off location of grids on next finer level
      !        c) Integrate over the grid and add to this procs total
      !     2) Accumulate/sum integral values from all procs
      !
      !
      locint  = 0.d0
      numpts  = 0
      lev     = level
  10  if (level_return_existence(lev)) then
         gi = level_return_start(lev)
  15     if (grid_return_existence(gi)) then
            if (grid_is_local(gi)) then
               if(ltrace)write(*,99)myid,'call grid_intfind: ',gi
               locint = locint+grid_intfind(gi,field,intorl2norm,tmppts)
               numpts = numpts + tmppts
            end if
            gi = grid_return_sibling(gi)
            goto 15
         end if
         if(ltrace)write(*,99)myid,'Done w/ lev: ',lev,locint
         lev = lev + 1
         goto 10
      end if
      !
      if (ltrace) then
         write(*,99)myid,'level_comp_integral: Local integral found'
         write(*,98)myid,'level_comp_integral: locint = ',locint
      end if
      !
      ! Get the global answer:
      call MPI_Reduce(locint,   globint, 1, MPI_DOUBLE_PRECISION,
     *                MPI_SUM,  master,  MPI_COMM_WORLD, ierr)
      if (intorl2norm .ne. INTEGRATE) then
         call MPI_Reduce(numpts,globnum, 1, MPI_INTEGER,
     *                MPI_SUM,  master,  MPI_COMM_WORLD, ierr)
         if(ltrace)write(*,99)myid,'level_comp_integralFinished',globnum
         if (globnum.ne.0) globint = sqrt(globint / globnum)
      end if
      if (ltrace)write(*,98)myid,'level_comp_integral:Finished',globint

      if (myid.eq.master) then
         call int2str(field,fieldnum)
         call grid_get_fname(field,fname)
         mylen           = mystringlength(fname)
         if (intorl2norm .eq. INTEGRATE) then
            fname(mylen+1:) = '_integral'
         else
            fname(mylen+1:) = '_l2norm'
         end if
         !write(*,*) 'Outputting integral to file: ',fname
         !open(unit=10,file=fname,STATUS='UNKNOWN')
         ! Only the master has the correct answer:
         !write(10,*) time, globint
#ifdef AIX
         gft_rc = gft_out_bbox( trim(fname)/ /CHAR(0), %VAL(time),
     *                    1, %VAL(1), bbox, globint )
#else
         gft_rc = gft_out_bbox( fname, time, 1, 1, bbox, globint )
#endif
      end if

  98  format('[',I3,'] ',A,F15.8)
  99  format('[',I3,'] ',A,3I5)
      return 
      end    ! END: level_comp_integral

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_comp_moments:                                                       cc
cc            Compute first and second moments, especially for                cc
cc            studying differentially rotating stars which may develop        cc
cc            a bar formation.                                                cc
cc         NB: routine parallels level_comp_integral since these moments      cc
cc             are just various integrals.                                    cc
cc         NB: See the following two references:                              cc
cc             http://arxiv.org/abs/astro-ph/0010201                          cc
cc             http://arxiv.org/abs/astro-ph/0609473                          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_comp_moments(level,field)
      implicit     none
      integer      level, field
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      lev, gi, i, j, mylen, gft_rc
      real(kind=8) time, bbox(2), tmp
      !
      ! Storing both second and first moments in same array:
      !       ( I^xx      I^xy     I^xz )
      !       ( I^x       I^yy     I^yz )
      !       ( I^y       I^z      I^zz )
      !
      real(kind=8)  moments(3,3)
      real(kind=8) gmoments(3,3)
      real(kind=8) etaplus, etacross
      data         bbox / 0.d0, 0.d0 /
      character(128) fname,  procname, fieldnum
      integer      gft_out_bbox, mystringlength
      external     gft_out_bbox, mystringlength

      logical      ltrace
      parameter (  ltrace = .false. )

      if (ltrace)write(*,99) myid, 'level_comp_moments: level: ',level

      time = grid_return_time(level_return_start(level))
      !
      !  We need to:
      !     1) For each local grid on this level:
      !        a) mask off location of other grids on level
      !        b) mask off location of grids on next finer level
      !        c) Integrate over the grid and add to this procs total
      !     2) Accumulate/sum integral values from all procs
      !
      !
      do j = 1, 3
         do i = 1, 3
         moments(i,j) = 0.d0
      end do
      end do
      lev     = level
  10  if (level_return_existence(lev)) then
         gi = level_return_start(lev)
  15     if (grid_return_existence(gi)) then
            if (grid_is_local(gi)) then
               if(ltrace)write(*,99)myid,'call grid_momentfind: ',gi
               call grid_momentfind( moments, gi, field)
            end if
            gi = grid_return_sibling(gi)
            goto 15
         end if
         if(ltrace)write(*,99)myid,'Done w/ lev: ',lev
         lev = lev + 1
         goto 10
      end if
      !
      if (ltrace) then
         write(*,99)myid,'level_comp_moments: Local integrals found'
         write(*,98)myid,'level_comp_moments: moments(1):',
     *        moments(1,1),moments(1,2),moments(1,3)
         write(*,98)myid,'level_comp_moments: moments(2):',
     *        moments(2,1),moments(2,2),moments(2,3)
         write(*,98)myid,'level_comp_moments: moments(3):',
     *        moments(3,1),moments(3,2),moments(3,3)
      end if
      !
      ! Get the global answer:
      call MPI_Reduce(moments,gmoments, 9, MPI_DOUBLE_PRECISION,
     *                MPI_SUM,  master,  MPI_COMM_WORLD, ierr)
      if (ltrace) then
         write(*,98)myid,'level_comp_moments:Finished:'
         if (myid.eq.master) then
            write(*,98)myid,'level_comp_moments:gmoments(1):',
     *           gmoments(1,1),gmoments(1,2),gmoments(1,3)
            write(*,98)myid,'level_comp_moments: moments(2):',
     *           gmoments(2,1),gmoments(2,2),gmoments(2,3)
            write(*,98)myid,'level_comp_moments:gmoments(3):',
     *           gmoments(3,1),gmoments(3,2),gmoments(3,3)
         end if
      end if

      ! compute distortion parameters:
      etaplus = (gmoments(1,1)-gmoments(2,2))
     *           /(gmoments(1,1)+gmoments(2,2))
      etacross = 2.d0*gmoments(1,2) /(gmoments(1,1)+gmoments(2,2))

      if (ltrace.and.myid.eq.master) then
         write(*,98)myid,'level_comp_moments:eta/plus/cross:',
     *           sqrt(etaplus**2+etacross**2),etaplus,etacross
      end if

      if (myid.eq.master) then
         call int2str(field,fieldnum)
         call grid_get_fname(field,fname)
         mylen           = mystringlength(fname)
         do j = 1,4
            do i = 1,3
               if (j.lt.4) tmp = gmoments(i,j)
               ! Output second moments:
               if      (i.eq.1 .and. j.eq.1) then
                  fname(mylen+1:) = '_xxmom'
               else if (i.eq.1 .and. j.eq.2) then
                  fname(mylen+1:) = '_xymom'
               else if (i.eq.1 .and. j.eq.3) then
                  fname(mylen+1:) = '_xzmom'
               else if (i.eq.2 .and. j.eq.2) then
                  fname(mylen+1:) = '_yymom'
               else if (i.eq.2 .and. j.eq.3) then
                  fname(mylen+1:) = '_yzmom'
               else if (i.eq.3 .and. j.eq.3) then
                  fname(mylen+1:) = '_zzmom'
               ! Output first moments:
               else if (i.eq.2 .and. j.eq.1) then
                  fname(mylen+1:) = '_xmom'
               else if (i.eq.3 .and. j.eq.1) then
                  fname(mylen+1:) = '_ymom'
               else if (i.eq.3 .and. j.eq.2) then
                  fname(mylen+1:) = '_zmom'
               ! Output distortion parameters:
               else if (i.eq.1 .and. j.eq.4) then
                  fname(mylen+1:) = '_etaplus'
                  tmp             = etaplus
               else if (i.eq.2 .and. j.eq.4) then
                  fname(mylen+1:) = '_etacross'
                  tmp             = etacross
               else if (i.eq.3 .and. j.eq.4) then
                  fname(mylen+1:) = '_eta'
                  tmp             = sqrt(etaplus**2+etacross**2)
               end if
#ifdef AIX
         gft_rc = gft_out_bbox( trim(fname)/ /CHAR(0), %VAL(time),
     *                    1, %VAL(1), bbox, tmp          )
#else
         gft_rc = gft_out_bbox( fname, time, 1, 1, bbox, tmp)
#endif
            end do
         end do
      end if

  98  format('[',I3,'] ',A,3F15.8)
  99  format('[',I3,'] ',A,3I5)
      return 
      end    ! END: level_comp_moments

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_modifydata_local:                                                   cc
cc            Allow the user to modify data, say after reading from           cc
cc            checkpoint, by calling a routine for each project found         cc
cc            in the comp_amr_error.f file.                                   cc
cc                                                                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_modifydata_local()
      implicit none
      include 'grid_methods.inc'
      integer  proc_return_myid
      external proc_return_myid
      integer  gi,     num,   lev_i, myid, nx,ny,nz, numspaces
      real(kind=8)      levwork,levh

      logical        ltrace
      parameter    ( ltrace = .false. )


      myid  = proc_return_myid()
      if(ltrace)write(*,99)myid,'level_modifydata_local: Enter.'
      if(ltrace.and.myid.eq.0) call level_tree_dump()
      lev_i = 0
      num   = 0
      if(ltrace)write(*,99)myid,'Consider      level: ',lev_i
  2   if (level_return_existence(lev_i)) then
         if(ltrace)write(*,99)myid,'Level exists:        ',lev_i
         gi        = level_return_start(lev_i)
         if(ltrace)write(*,99)myid,'Level starts with gi:',gi
  5      if (grid_return_existence(gi)) then
            if(ltrace)write(*,99)myid,'Grid exists:         ',gi
            if ( grid_is_local(gi) ) then
               if(ltrace)write(*,99)myid,'Grid local:          ',gi
               call grid_modify_data(gi,nx,ny,nz)
            end if
            gi = grid_return_sibling(gi)
            if(ltrace)write(*,99)myid,'Consider next grid:  ',gi
            goto 5
         end if
         lev_i = lev_i + 1
         if(ltrace)write(*,99)myid,'Consider next level: ',lev_i
         num   = 0
         goto 2
      end if
      if(ltrace) call level_test_local()
      if(ltrace) write(*,99)myid,'level_modifydata_local: Done.'

  98  format('[',I4,'] ',A,3F15.8)
  99  format('[',I4,'] ',A,3I5)
      return
      end      ! END: subroutine level_modifydata_local


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_test_local:                                                         cc
cc            For debugging...call grid_test() on each local grid.            cc
cc                                                                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_test_local()
      implicit none
      include 'grid_methods.inc'
      integer  proc_return_myid
      external proc_return_myid
      integer  gi,     num,   lev_i, myid, nx,ny,nz, numspaces
      real(kind=8)      levwork,levh

      logical        ltrace
      parameter    ( ltrace = .false. )


      myid  = proc_return_myid()
      if(ltrace)write(*,99)myid,'level_test_local: Enter.'
      lev_i = 0
      num   = 0
      if(ltrace)write(*,99)myid,'Consider      level: ',lev_i
  2   if (level_return_existence(lev_i)) then
         if(ltrace)write(*,99)myid,'Level exists:        ',lev_i
         gi        = level_return_start(lev_i)
         if(ltrace)write(*,99)myid,'Level starts with gi:',gi
  5      if (grid_return_existence(gi)) then
            if(ltrace)write(*,99)myid,'Grid exists:         ',gi
            if ( grid_is_local(gi) ) then
               if(ltrace)write(*,99)myid,'level_test:Grid local ',gi
               call grid_test(gi)
            end if
            gi = grid_return_sibling(gi)
            if(ltrace)write(*,99)myid,'level_test:Consider next grid',gi
            goto 5
         end if
         lev_i = lev_i + 1
         if(ltrace)write(*,99)myid,'level_test:Consider nextlevel',lev_i
         num   = 0
         goto 2
      end if
      if(ltrace)write(*,99)myid,'level_test_local: Done.'

  98  format('[',I4,'] ',A,3F15.8)
  99  format('[',I4,'] ',A,3I5)
      return
      end      ! END: subroutine level_test_local

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  level_central_value:                                                      cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine level_central_value(level,field)
      implicit     none
      integer      level, field
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      MIN,     MAX,     ABSMIN,     ABSMAX
      parameter (  MIN = 0, MAX = 1, ABSMIN = 2, ABSMAX = 3 )
      integer      lev, gi, i
      integer      mpi_opt
      ! Processor ID of who owns the min or max:
      integer      procext
      real(kind=8) locext,  loccoords(3)
      real(kind=8) globext, globcoords(3)
      real(kind=8)   LARGENUMBER
      parameter    ( LARGENUMBER = 9.d98 )
      logical      double_equal
      external     double_equal
      ! Output variables:
      integer      mylen, gft_rc, rank, shape(3)
      parameter  ( rank = 3)
      data         shape / 1, 1, 1 /
      real(kind=8) time, bbox(6)
      character(128) fname,  procname, fieldnum
      ! Count the number of points over which l2norm is taken
      integer      tmppts, numpts, globnum
      integer      gft_out_bbox, mystringlength
      external     gft_out_bbox, mystringlength

      logical      ltrace
      parameter (  ltrace  = .false. )
      logical      ltrace2
      parameter (  ltrace2 = .false. )

      time = grid_return_time(level_return_start(level))
      if (ltrace) then
         write(*,99) myid, '----------------------------------'
         write(*,99) myid, 'level_central_value: level: ',level
      end if

      !
      lev     = level
  10  if (level_return_existence(lev)) then
         gi = level_return_start(lev)
  15     if (grid_return_existence(gi)) then
            if (grid_is_local(gi)) then
               if(ltrace2)write(*,99)myid,'call grid_extfind: ',gi
               !call grid_extfind(gi, field, locext,loccoords, minormax)
            end if
            gi = grid_return_sibling(gi)
            goto 15
         end if
         lev = lev + 1
         goto 10
      end if
      !
      if (ltrace) then
         write(*,99)myid,'level_comp_extremum: Local min/max found'
         write(*,98)myid,'level_comp_extremum: locext:',locext
         do i = 1, 3
            write(*,98)myid,'level_comp_extremum: coords: ',loccoords(i)
         end do
      end if
      !

      ! Tell all procs the answer:
      call MPI_BCAST(globext, 1, MPI_DOUBLE_PRECISION,
     &               master, MPI_COMM_WORLD, ierr)
      if(ltrace)write(*,98)myid,' globext: ',globext


      do i = 1, 3
         bbox(2*i-1) = loccoords(i)
         bbox(2*i  ) = loccoords(i)
      end do
      !
      ! Output the min or max:
      !
      if (myid.eq.master) then
         call int2str(field,fieldnum)
         call grid_get_fname(field,fname)
         mylen           = mystringlength(fname)
         fname(mylen+1:) = '_central'
         if (ltrace) then
            write(*,99)myid,' procext:  ',procext
            write(*,99)myid,' numprocs: ',numprocs
            write(*,98)myid,' locext:   ',locext
            write(*,98)myid,' globext:  ',globext
            write(*,97)myid,' Coords:   ',1,bbox(1)
            write(*,97)myid,' Coords:   ',2,bbox(3)
            write(*,97)myid,' Coords:   ',3,bbox(5)
            write(*,*) myid,' fname: ',fname(1:13)
         end if
#ifdef AIX
         gft_rc = gft_out_bbox( trim(fname)/ /CHAR(0), %VAL(time),
     *                    shape, %VAL(rank), bbox, globext )
#else
         gft_rc = gft_out_bbox( fname, time, shape, rank, bbox, globext)
#endif
      end if

      if (ltrace) then
         write(*,98) myid, 'level_comp_extremum: Finished',globext
         do i = 1, 3
            write(*,97) myid, 'level_comp_extremum: ',i,loccoords(i)
         end do
         write(*,99) myid, '----------------------------------'
      end if

  97  format('[',I3,'] ',A,I4,F15.8)
  98  format('[',I3,'] ',A,F15.8)
  99  format('[',I3,'] ',A,3I5)
      return 
      end    ! END: level_comp_extremum
