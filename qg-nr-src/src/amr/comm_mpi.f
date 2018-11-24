cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  proc_pick_group:                                                          cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine proc_pick_group(bminx,bmaxx,bminy,bmaxy,bminz,bmaxz,
     *                           childowners, minxg, minyg, minzg, 
     *                           hi, nbox, level)
      implicit     none
      integer      nbox, level
      integer      bminx(nbox), bmaxx(nbox),
     *             bminy(nbox), bmaxy(nbox),
     *             bminz(nbox), bmaxz(nbox),
     *             childowners(nbox)
      real(kind=8) minxg,minyg,minzg, hi
      include     'glob.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      real(kind=8) minx, miny, minz, maxx, maxy, maxz
      integer      i, j, size, numgfuncs, gi, proc
      integer      spaceavailable, orig_pick,spaceneed
      integer      proc_overlap(0:maxprocs-1)
      integer      proc_olap_max, owner, tmpi
      integer      proc_pick
      !
      integer      grid_bbox_intersection
      external     grid_bbox_intersection
      !
      integer      proc_least_loaded
      external     proc_least_loaded
      !
      integer      mem_return_arena
      external     mem_return_arena
      !
      logical      ltrace
      parameter  ( ltrace = .true. )


      if (ltrace) then
         write(*,*) 'proc_pick_group: Choosing owners for nbox: ',nbox
         write(*,*) 'proc_pick_group: On level:                 ',level
         call proc_dump_wkld()
      end if

      if (level.lt.0) then
         if (ltrace) write(*,*) 'proc_pick_group: Going sequentially'
         !
         ! Creating coarse grid(s):
         !    (don't need to worry about issues)
         !
         proc = 0
         do i = 1, nbox
            childowners(i) = proc
            proc           = proc + 1
            if (proc.eq.numprocs) proc = 0
         end do
         return
      end if

      numgfuncs = grid_return_numgfuncs()

      !
      ! Loop over the different bounding boxes and assign owners:
      !
      do i = 1, nbox
         !
         ! Compute space needed:
         size           =   ((bmaxx(i)-bminx(i))*refine_factor+1)
     *                     *((bmaxy(i)-bminy(i))*refine_factor+1)
     *                     *((bmaxz(i)-bminz(i))*refine_factor+1)
         spaceneed      = size * numgfuncs
         !
         minx = minxg + (bminx(i)-1) * hi
         miny = minyg + (bminy(i)-1) * hi
         minz = minzg + (bminz(i)-1) * hi
         maxx = minxg + (bmaxx(i)-1) * hi
         maxy = minyg + (bmaxy(i)-1) * hi
         maxz = minzg + (bmaxz(i)-1) * hi
         !
         ! Determine which proc owns grids with most overlap:
         !
         ! Initialize proc info:
         proc_olap_max = -1
         do j = 0, maxprocs-1
            proc_overlap(j) = 0
         end do
         ! Loop over level to be refined to determine overlaps:
         if (ltrace) write(*,*)'proc_pick_group: Searching most overlap'
         gi = level_return_start(level)
 10      if (grid_return_existence(gi)) then
            owner = grid_return_owner(gi)
            tmpi  = grid_bbox_intersection(gi,minx,maxx,miny,
     *                                        maxy,minz,maxz)
            proc_overlap(owner) = proc_overlap(owner) + tmpi
            if (proc_olap_max .lt. 0 .or.
     *          proc_overlap(owner) .gt. proc_overlap(proc_olap_max) )
     *          proc_olap_max = owner
            if (ltrace) then
               write(*,*)'proc_pick_group: Examining grid: ',gi
               write(*,*)'proc_pick_group: with overlap:   ',tmpi
               write(*,*)'proc_pick_group: with owner:     ',owner
            end if
            gi = grid_return_sibling(gi)
            goto 10
         end if
         proc_pick = proc_olap_max
         !
         ! Determine whether to overule this choice:
         !
         ! check if enough memory available on this proc:
         spaceavailable = mem_return_arena() -NINT(proc_meml(proc_pick))
         if (spaceneed .ge. spaceavailable) then
            write(*,*)'proc_pick_group: Not enough space: need/avail ',
     *                                          spaceneed,spaceavailable
         end if
         childowners(i) = proc_pick
         if (ltrace) then
            write(*,*) 'proc_pick_group: Owner for box ',i,
     *                                      ' is ',childowners(i)
            write(*,*) 'proc_pick_group: Size of box i: ',size
         end if
      end do


      if (ltrace) then
         write(*,*) 'proc_pick_group: Finished choosing.'
      end if

      return
      end          ! END: proc_pick_group

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  proc_pick:                                                                cc
cc             First try at finding a good processor on which                 cc
cc             to create a grid. Called by grid_refine().                     cc
cc             Takes as input the grid that is getting refined.               cc
cc             Updated: to pick minimum with sufficient memory.               cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function proc_pick(grid, size)
      implicit     none
      integer      grid, size
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'grid_methods.inc'
      integer      spaceavailable, orig_pick,spaceneed
      integer      totalarena
      ! Rough estimate of overhead memory (not associated with fields):
      !                                   (e.g. x/y/z arrays)            
      integer      OVERHEADMEM
      parameter (  OVERHEADMEM = 500 000 )
      !
      logical      nogridsonmaster
      parameter  ( nogridsonmaster  = .false. )
      !
      integer      proc_least_loaded
      external     proc_least_loaded
      !
      logical      proc_havemem
      external     proc_havemem
      !
      integer      mem_return_arena, mem_return_avail
      external     mem_return_arena, mem_return_avail
      !
      integer      find_min_index1D
      external     find_min_index1D

      logical      ltrace
      parameter  ( ltrace = .false. )

 97   format('[',I4,']',A,E16.7)
 98   format('[',I4,']',A)
 99   format('[',I4,']',A,3I12)

      !if (ltrace) call proc_dump_wkld()

      !
      ! Size is size of grid, to get memory needs,
      ! need to multiply by the number of grid functions
      !
      spaceneed = size * grid_return_numgfuncs()

      if (ltrace) write(*,99) myid,'proc_pick: Pick proc for grid ',grid
      if (ltrace) write(*,99) myid,'proc_pick: of size:           ',size
      if (ltrace) write(*,99) myid,'proc_pick: with space need of ',
     *      spaceneed

      proc_pick = proc_least_loaded()
      if (ltrace) write(*,99) myid,'proc_pick: least loaded: ',proc_pick

      if(nogridsonmaster.and.proc_pick.eq.0 .and.numprocs.gt.1) then
         proc_pick = 1
         if (ltrace) write(*,99) myid,'proc_pick: do not go with 0'
      end if

      !
      ! Does processor have enough memory for the grid?
      !  (assuming all procs have the same arena size)
      !
      orig_pick      = proc_pick

c     totalarena     = mem_return_arena()
c     spaceavailable = totalarena - NINT(proc_meml(proc_pick))
c    *                            - OVERHEADMEM
c     ! If local, we can get precise amount of memory available:
c     if (myid .eq. proc_pick) spaceavailable = mem_return_avail()
c    *                            - OVERHEADMEM
c     if (ltrace) write(*,99) myid,'proc_pick: space available:   ',
c    *                                                    spaceavailable

      if (ltrace) call mem_stat()
      if (.not. proc_havemem(proc_pick,spaceneed) ) then
      !if (spaceneed .ge. spaceavailable) then
         write(*,*)'proc_pick: Cannot place on processor: ',proc_pick
         write(*,*)'proc_pick: Not enough space: need/avail ',
c    *                                          spaceneed,spaceavailable
     *                                          spaceneed
         ! Because of indexing from 0
         proc_pick = find_min_index1D(proc_meml(0),numprocs) - 1
         write(*,*)'proc_pick: Trying the least memloaded proc: ',
     *                                          proc_pick
c        spaceavailable = totalarena - NINT(proc_meml(proc_pick))
c    *                               - OVERHEADMEM
c        if (myid .eq. proc_pick) spaceavailable = mem_return_avail()
c    *                               - OVERHEADMEM
         if (.not. proc_havemem(proc_pick,spaceneed) ) then
         !if (spaceneed .gt. spaceavailable) then
            write(*,*)'proc_pick: Leastloaded cannot: ',proc_pick
            write(*,*)'proc_pick: So no processor available ****.'
            write(*,*)'proc_pick: Not enough space: need/avail ',
c    *                                          spaceneed,spaceavailable
     *                                          spaceneed
            proc_pick = -1
            call proc_dump_wkld()
            return
         end if
      end if
      if (proc_pick .ne. orig_pick) then
         write(*,*) 'proc_pick: Done.  ------Chose: ',proc_pick
         call proc_dump_wkld()
      end if

      if (ltrace) write(*,99) myid,'proc_pick: Done, choose ',proc_pick

      return
      end          ! END: proc_pick

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  proc_havemem:                                                             cc
cc                   Returns true if processor has enough memory else false.  cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function proc_havemem(proc,size)
      implicit     none
      integer      proc, size
      include     'mpif.h'
      include     'mpi_stuff.inc'
      integer      totalarena, spaceavailable
      !
      integer      OVERHEADMEM
      parameter (  OVERHEADMEM = 1 000 000 )
      !
      integer      mem_return_arena, mem_return_avail
      external     mem_return_arena, mem_return_avail

      totalarena     = mem_return_arena()

      spaceavailable = totalarena - NINT(proc_meml(proc))
     *                            - OVERHEADMEM
      ! If local, we can get precise amount of memory available:
      if (myid .eq. proc) spaceavailable = mem_return_avail()
     *                            - OVERHEADMEM

      if (size .ge. spaceavailable) then
         proc_havemem = .false.
      else
         proc_havemem = .true.
      end if

      return
      end          ! END: proc_havemem


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  proc_least_loaded:                                                        cc
cc                     Return processor id with least load.                   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function proc_least_loaded()
      implicit     none
      include     'mpif.h'
      include     'mpi_stuff.inc'

      proc_least_loaded = proc_minwld_id

      return
      end          ! END: proc_least_loaded

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  barrier:                                                                  cc
cc              Called by master process, to ensure that                      cc
cc              other processes are caught up with commands.                  cc
cc              It suffices to send each processor an instruction             cc
cc              to respond, since this instruction will come after            cc
cc              the others from the master process.                           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine barrier()
      implicit     none
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'action.inc'
      integer      dest, tmp_int
      integer      request(1:numprocs-1)

      logical      ltrace
      parameter  ( ltrace = .false. )

    
      if (ltrace) write(*,99) myid, 'barrier: Enter'

      if (myid .ne. master) then
         write(*,99) myid, 'barrier: Must be master to run'
         call my_exit('Error in my_barrier')
      end if

      !
      ! Request acknowledgement from other procs:
      !
      if (ltrace) write(*,99) myid, 'barrier: Sending MY_BARRIER'
      call send_actionAll(MY_BARRIER, 0, 0,0,0, 0,0,0, 0,0,0)

      do dest = 1, numprocs - 1
         !
         ! Receive acknowledgement from other procs:
         !
         if (ltrace) write(*,99) myid, 'barrier: IRecv for: ',dest
         call MPI_IRecv(tmp_int, 1, MPI_INTEGER,
     *           dest, TAG_BARRIER, MPI_COMM_WORLD,request(dest),ierr)
      end do

      do dest = 1, numprocs - 1
         if (ltrace) write(*,99) myid, 'barrier: Wait  for: ',dest
         call MPI_Wait(request(dest), status, ierr)
      end do

      if (ltrace) write(*,99) myid, 'barrier: Exit'

 99   format('[',I4,'] ',A,I3)

      return
      end          ! END: barrier

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  send_action:                                                              cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine send_action(dest, action, gi, owner, level, parent,
     *                        gnx,gny,gnz,gax,gay,gaz)
      implicit     none
      integer      dest, action, gi, owner, level, parent,
     *             gnx, gny, gnz, gax, gay, gaz
      include     'mpif.h'
      include     'mpi_stuff.inc'

      logical      ltrace
      parameter  ( ltrace = .false. )

    
      !
      !  Check destination processor:
      !
      if (dest.lt.0 .or. dest.ge.numprocs) return

      if(ltrace)write(*,91)myid,'send_action: action = ',action

         gridmsg(1)  = action
         gridmsg(2)  = gi
         gridmsg(3)  = owner
         gridmsg(4)  = level
         gridmsg(5)  = parent
         gridmsg(6)  = gnx
         gridmsg(7)  = gny
         gridmsg(8)  = gnz
         gridmsg(9)  = gax
         gridmsg(10) = gay
         gridmsg(11) = gaz

      if(ltrace)write(*,91)myid,'send_action: sending to: ',dest
      call MPI_Send(gridmsg, mlength, MPI_INTEGER,
     *                 dest, TAG_ACTION, MPI_COMM_WORLD, ierr)

      if(ltrace)write(*,90)myid,'send_action: Action sent.'

 90   format('[',I4,'] ',A)
 91   format('[',I4,'] ',A,I3)
 99   format('[',I4,'] ',A,I3)

      return
      end          ! END: send_action

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  send_actionOthers:                                                        cc
cc                 Brute force to send to all other processors.               cc
cc              NB: does *not* send to the master process.                    cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine send_actionOthers(action, gi, owner, level, parent,
     *                        gnx,gny,gnz,gax,gay,gaz)
      implicit     none
      integer      action, gi, owner, level, parent,
     *             gnx, gny, gnz, gax, gay, gaz
      include     'mpif.h'
      include     'mpi_stuff.inc'
      integer      dest, request(1:numprocs-1)

      logical      ltrace
      parameter  ( ltrace  = .true. )

    
         gridmsg(1)  = action
         gridmsg(2)  = gi
         gridmsg(3)  = owner
         gridmsg(4)  = level
         gridmsg(5)  = parent
         gridmsg(6)  = gnx
         gridmsg(7)  = gny
         gridmsg(8)  = gnz
         gridmsg(9)  = gax
         gridmsg(10) = gay
         gridmsg(11) = gaz

      do dest = 1, numprocs - 1
         if (dest .ne. myid) then
           if (ltrace) write(*,97) myid,'send_actionOthers: dest:',dest
           call MPI_ISend(gridmsg, mlength, MPI_INTEGER,
     *            dest, TAG_ACTION, MPI_COMM_WORLD, request(dest),ierr)
         end if
      end do

      do dest = 1, numprocs - 1
         if (dest .ne. myid) then
            call MPI_Wait(request(dest), status, ierr)
         end if
      end do

      if (ltrace) write(*,99) myid,'send_actionOthers: Action sent.'

 97   format('[',I4,'] ',A,2I3)
 98   format('[',I4,'] ',A, I3)
 99   format('[',I4,'] ',A)

      return
      end          ! END: send_actionOthers

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  send_actionAllslow:                                                       cc
cc                 loop over all processes and send an action to each.        cc
cc                 This is called 'slow' in comparison to the non-slow        cc
cc                 one because that one had problems with incorrect           cc
cc                 sequencing of commands.                                    cc
cc              NB: can be called by any processor.                           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine send_actionAllslow(action, gi, owner, level, parent,
     *                        gnx,gny,gnz,gax,gay,gaz)
      implicit     none
      integer      action, gi, owner, level, parent,
     *             gnx, gny, gnz, gax, gay, gaz
      include     'mpif.h'
      include     'mpi_stuff.inc'
      integer      dest, request(0:numprocs)

      logical      ltrace
      parameter  ( ltrace  = .false. )

    
         gridmsg(1)  = action
         gridmsg(2)  = gi
         gridmsg(3)  = owner
         gridmsg(4)  = level
         gridmsg(5)  = parent
         gridmsg(6)  = gnx
         gridmsg(7)  = gny
         gridmsg(8)  = gnz
         gridmsg(9)  = gax
         gridmsg(10) = gay
         gridmsg(11) = gaz

      do dest = 0, numprocs - 1
         if (dest .ne. myid)
     *   call MPI_ISend(gridmsg, mlength, MPI_INTEGER,
     *            dest, TAG_ACTION, MPI_COMM_WORLD, request(dest),ierr)
         if (ltrace) write(*,97) myid,'           dest: ',dest
      end do

      do dest = 0, numprocs - 1
         if (dest .ne. myid)
     *   call MPI_Wait(request(dest), status, ierr)
      end do

      if (ltrace) write(*,99) myid,'send_actionAllslow: Action sent.'

 97   format('[',I4,'] ',A,2I3)
 98   format('[',I4,'] ',A, I3)
 99   format('[',I4,'] ',A)

      return
      end          ! END: send_actionAllslow

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  send_actionAll:                                                           cc
cc                 send a command to all other processors                     cc
cc                 To speed things up, processors which already               cc
cc                 have received this message are also used to send           cc
cc                 to some of the processors. Hence, one can think            cc
cc                 of the process as going in rounds. In the zeroth           cc
cc                 round, the master process (0) sends to processor 1.        cc
cc                 In the next round, 0 sends to 2 and 1 sends to 3,          cc
cc                 and so on:                                                 cc
cc                       0  0  0                                              cc
cc                          1  1                                              cc
cc                             2                                              cc
cc                             3                                              cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine send_actionAll(action, gi, owner, level, parent,
     *                        gnx,gny,gnz,gax,gay,gaz)
      implicit     none
      integer      action, gi, owner, level, parent,
     *             gnx, gny, gnz, gax, gay, gaz
      include     'mpif.h'
      include     'mpi_stuff.inc'
      integer      dest, request(0:numprocs)
      integer      round, start_round, i
      integer      log_basen_INT
      external     log_basen_INT

      logical      SLOW  
      parameter  ( SLOW    = .false. )

      logical      ltrace
      parameter  ( ltrace  = .false. )
      logical      ltrace2
      parameter  ( ltrace2 = .false. )

      if (SLOW) then
         if (myid .eq. master) then
            call send_actionAllslow(action, gi, owner, level, parent,
     *                        gnx,gny,gnz,gax,gay,gaz)
         end if
         if (ltrace) write(*,98) myid,'send_actionAll: Done w/ slow',
     *                    action
         return
      end if
    
      if(ltrace)write(*,98)myid,'send_actionAll: Sending action:',
     *                    action
         gridmsg(1)  = action
         gridmsg(2)  = gi
         gridmsg(3)  = owner
         gridmsg(4)  = level
         gridmsg(5)  = parent
         gridmsg(6)  = gnx
         gridmsg(7)  = gny
         gridmsg(8)  = gnz
         gridmsg(9)  = gax
         gridmsg(10) = gay
         gridmsg(11) = gaz

      !
      ! Determine in which round we enter this process:
      !
      if (myid .eq. 0) then
         start_round = 0
      else 
         start_round = 1 + log_baseN_INT(myid,2)
      end if
      round = start_round
      if (ltrace2) then
         write(*,98) myid, 'send_actionAll: Start_round: ',start_round
      end if

 10   dest = myid + 2**round
      if (dest .lt. numprocs) then
         call MPI_ISend(gridmsg, mlength, MPI_INTEGER,
     *            dest, TAG_ACTION, MPI_COMM_WORLD, request(round),ierr)
         if (ltrace) write(*,97) myid,'     round/dest: ',round, dest
         round = round + 1
         goto 10
      end if

      do i = start_round, round-1
         if (ltrace) write(*,97) myid,'  Waiting for round:',i
         call MPI_Wait(request(i), status, ierr)
      end do

      if (ltrace2) write(*,99) myid,'send_actionAll: Action sent.'

 97   format('[',I4,'] ',A,2I5)
 98   format('[',I4,'] ',A, I3)
 99   format('[',I4,'] ',A)

      return
      end          ! END: send_actionAll

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  recv_action:                                                              cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine recv_action(action, gi, owner, level, parent,
     *                        gnx,gny,gnz,gax,gay,gaz,from)
      implicit     none
      integer      action, gi, owner, level, parent,
     *             gnx, gny, gnz, gax, gay, gaz, from
      include     'mpif.h'
      include     'mpi_stuff.inc'
      integer      recv_error, recv_tag

      logical      ltrace
      parameter  ( ltrace = .false. )

      if(ltrace)write(*,99)myid,'recv_action: Receiving action'

         call MPI_RECV(gridmsg, mlength, MPI_INTEGER,
     *          MPI_ANY_SOURCE, TAG_ACTION, MPI_COMM_WORLD,status, ierr)

      !
      !  status(1)   ----SGI compilers and arrayd
      !  status(2)   ----On Linux with Intel Compilers & MPICH
      !              ----On Macs, with g77 &  MPICH
      !              ----On titan, with Intel 7.0 & MPICH
      !
      from       = status(MPI_SOURCE)
      recv_error = status(MPI_ERROR)
      recv_tag   = status(MPI_TAG)

      from = status(MPI_SOURCE)
c     from = status(2)
c     from = status(1)

      if (ltrace) then
         write(*,98)myid,'  ierr:         ',ierr
         write(*,98)myid,'  Message from: ',from
         write(*,98)myid,'  tag:          ',recv_tag
      end if

      !
      ! There seems to be some disparity among implemnations
      ! of MPI. The sender (when receiving with MPI_ANY_SOURCE)
      ! is supposed to be returned in status(). In some
      ! implementation it was in status(1) whereas now it's in 
      ! status(2). Keepin mind status(1..4).
      !
      !
      if(ltrace) then
         write(*,98) myid, '  Message from 1: ',status(1)
         write(*,98) myid, '  Message from 2: ',status(2)
         write(*,98) myid, '  Message from 3: ',status(3)
         write(*,98) myid, '  Message from 4: ',status(4)
      end if

      if (.false.) then
         write(*,*) gridmsg(1)
         write(*,*) gridmsg(2)
         write(*,*) gridmsg(3)
         write(*,*) gridmsg(4)
         write(*,*) gridmsg(5)
         write(*,*) gridmsg(6)
         write(*,*) gridmsg(7)
         write(*,*) gridmsg(8)
         write(*,*) gridmsg(9)
         write(*,*) gridmsg(10)
         write(*,*) gridmsg(11)
      end if

         action = gridmsg(1)
         gi     = gridmsg(2)
         owner  = gridmsg(3)
         level  = gridmsg(4)
         parent = gridmsg(5)
         gnx    = gridmsg(6)
         gny    = gridmsg(7)
         gnz    = gridmsg(8)
         gax    = gridmsg(9)
         gay    = gridmsg(10)
         gaz    = gridmsg(11)

      if(ltrace)write(*,98)myid,'recv_action: Action received.',action


 97   format('[',I4,'] ',A,2I3)
 98   format('[',I4,'] ',A, I3)
 99   format('[',I4,'] ',A)

      return
      end          ! END: recv_action


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  proc_subtractwork:                                                        cc
cc                 Subtract the work of a grid from  entry of processor       cc
cc                 in table. At same time, update which processor has         cc
cc                 the minimum workload and what that workload is.            cc
cc                 Updated: to remove memory load for the proc.               cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine proc_subtractwork(id, work, size)
      implicit     none
      integer      id, size
      real(kind=8)     work
      include     'mpif.h'
      include     'mpi_stuff.inc'
      integer      nextproc
      real(kind=8)     tmp_load, tmp_mload
      real(kind=8)     SMALLNUMBER
      parameter      ( SMALLNUMBER = 1.0 d-10)
      !
      integer      grid_return_numgfuncs
      external     grid_return_numgfuncs

      logical      nogridsonmaster
      parameter  ( nogridsonmaster  = .false. )
      logical      ltrace
      parameter  ( ltrace  = .false. )
      logical      ltrace2
      parameter  ( ltrace2 = .false. )

 97      format('[',I4,'] ',A,I4,3F16.8)
 98      format('[',I4,'] ',A,3G16.3)
 99      format('[',I4,'] ',A,3I15)
      if(ltrace2.and.myid.eq.master .or. ltrace2)
     *     write(*,97)myid,'proc_subtractwork: id, work ', id, work

      tmp_mload     = (grid_return_numgfuncs()+1)*size
      tmp_load      = proc_wkld(id) - work
      proc_wkld(id) = tmp_load
      if (tmp_load .lt. 0) then
         if (tmp_load.gt.SMALLNUMBER) then
            write(*,99)myid,' proc_subtractwork: Work cant be negative'
            write(*,99)myid,' proc_subtractwork:     id = ', id
            write(*,98)myid,' proc_subtractwork:   work = ', work
            write(*,98)myid,' proc_subtractwork: resulti= ', tmp_load
            write(*,99)myid,' proc_subtractwork: Setting to zero'
         end if
         proc_wkld(id) = 0.d0
         tmp_load      = 0.d0
      end if
      !proc_meml(id) = proc_meml(id) - tmp_mload
      if (proc_meml(id) .lt. tmp_mload) then
      !if (proc_meml(id) .lt. 0) then
          write(*,99)myid,' proc_subtractwork: Memload cant go negative'
          write(*,99)myid,' proc_subtractwork:     id = ', id
          write(*,98)myid,' proc_subtractwork:   work = ', work
          write(*,99)myid,' proc_subtractwork:   size = ', size
          write(*,98)myid,' proc_subtractwork: proc_meml(id) = ',
     *                                         proc_meml(id)
          write(*,98)myid,' proc_subtractwork: tmp_mload = ',tmp_mload
          call proc_dump_wkld()
          proc_meml(id) = 0.d0
      else
          proc_meml(id) = proc_meml(id) - tmp_mload
      end if

      if (ltrace2) write(*,*) 'proc_subtractwork: Subtracted workload'
      if (ltrace2) call proc_dump_wkld()
      if (ltrace .and. myid .eq. master) then
         write(*,91)id,nint(tmp_mload),proc_wkld(id),NINT(proc_meml(id))
 91      format('-             id: ',I4,' mload: ',i13,' wkld: ',f9.3,
     *                       ' meml: ',I13)
      end if

      !
      ! Keep track of minimum loaded processor:
      !
      if ( .not.(nogridsonmaster .and. id.eq.0) ) then
         if (proc_minwld_id .ne. id .and. numprocs .ne. 1) then
            !
            ! See if this is now least loaded
            !
            if (tmp_load .lt. proc_minwld) then
               proc_minwld_id = id
               proc_minwld    = proc_wkld(id)
            end if
         else
            !
            ! This remains least loaded proc, but with less load:
            !
            proc_minwld    = proc_wkld(id)
         end if
      end if

      if (ltrace2)
     *   write(*,*) 'proc_subtractwork: minwld_id:',proc_minwld_id
      if (ltrace2)
     *   write(*,*) 'proc_subtractwork: minwld:   ',proc_minwld

 20   return
      end          ! END: proc_subtractwork

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  proc_addwork:                                                             cc
cc                 Add an amount of work to workload entry of processor       cc
cc                 in table. At same time, update which processor has         cc
cc                 the minimum workload and what that workload is.            cc
cc                 Updated: to include storage of memory load.                cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine proc_addwork(id, work, size)
      implicit     none
      integer      id, size
      real(kind=8)     work
      include     'mpif.h'
      include     'mpi_stuff.inc'
      integer      nextproc
      real(kind=8)     tmp_load, tmp_mload
      !
      logical      nogridsonmaster
      parameter  ( nogridsonmaster  = .false. )
      !
      integer      grid_return_numgfuncs
      external     grid_return_numgfuncs

      logical      ltrace
      parameter  ( ltrace  = .false. )
      logical      ltrace2
      parameter  ( ltrace2 = .false. )

      if(ltrace2.and.myid.eq.master .or. ltrace2)
     *    write(*,*)myid,'proc_addwork: id, work,size ', id, work, size

      !
      ! Keep track of memory load on each proc:
      !    (add one for storage of x/y/z arrays)
      !
      tmp_mload     = (grid_return_numgfuncs()+10)*size
      !tmp_mload     = (grid_return_numgfuncs()+1)*size
      tmp_load      = proc_wkld(id) + work
      proc_wkld(id) = tmp_load
      proc_meml(id) = proc_meml(id) + tmp_mload

      if(ltrace.and.myid.eq.master .or. ltrace2) then
          !write(*,*) 'proc_addwork: Added workload'
          !call proc_dump_wkld()
         write(*,91)myid,id,nint(tmp_mload),proc_wkld(id),
     &                NINT(proc_meml(id))
 91      format(I4']+             id: ',I4,' mload: ',i13,' wkld: ',f9.3
     *                       ,' meml: ',I13)
      end if

      !
      ! Keep track of minimum loaded processor:
      !
      !if (proc_minwld_id .eq. id .and. numprocs .ne. 1) then
      if ( (proc_minwld_id .eq. id .and. numprocs .ne. 1) 
     *                    .or.
     *      (proc_minwld_id .eq.0 .and. nogridsonmaster)   )  then
         if (id.eq.0.and.nogridsonmaster) then
            ! Do not allow master to be considered min loaded:
            ! and so:   (i)  set minloaded to proc 1
            !           (ii) then search if others less loaded exc. master
            proc_minwld_id = 1
            tmp_load       = proc_wkld(1)
            proc_minwld    = tmp_load
         end if
         !
         ! Find new minimum:
         !
         nextproc = id + 1
         if (nextproc .eq. numprocs) nextproc = 0
         if(nogridsonmaster.and.nextproc.eq.0) nextproc=1
 10      continue
         if (ltrace2)write(*,*)myid,'  nextproc:  ',nextproc
         if (ltrace2)write(*,*)myid,'  w/ load:   ',proc_wkld(nextproc)
         if (ltrace2)write(*,*)myid,'  numprocs:  ',numprocs
         if (ltrace2)write(*,*)myid,'  id:        ',id
         if ( proc_wkld(nextproc) .eq. 0 ) then
            !
            ! No need to go through whole list:
            !
            proc_minwld_id = nextproc
            proc_minwld    = proc_wkld(nextproc)
            goto 20
         else if ( proc_wkld(nextproc) .lt. tmp_load ) then
            proc_minwld_id = nextproc
            tmp_load       = proc_wkld(nextproc)
            proc_minwld    = tmp_load
         end if
         nextproc = nextproc + 1
         if (nextproc .eq. numprocs) nextproc = 0
         if(nogridsonmaster.and.nextproc.eq.0.and.id.ne.0) nextproc=1
         if (nextproc .ne. id) goto 10
      end if

      if(ltrace2.and.myid.eq.master .or. ltrace2) then
      !if(ltrace.and.myid.eq.master .or. ltrace2) then
         write(*,90) myid,proc_minwld_id, proc_minwld
         !write(*,*) 'proc_addwork: minwld_id:',proc_minwld_id
         !write(*,*) 'proc_addwork: minwld:   ',proc_minwld
 90      format(I4,']proc_addwork: minwld_id: ',I4,'  minwld: ',F9.3)
      end if

 20   return
      end          ! END: proc_addwork

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  proc_addmem:                                                              cc
cc                 Add an amount of memory used by a process to its entry.    cc
cc                 Note that in contrast to proc_addwork, this does *not*     cc
cc                 assume that this is the amount of memory associated        cc
cc                 with each field on a grid. Instead, it is just a set       cc
cc                 amount of memory (as would be used by the master process   cc
cc                 in the refinement process).                                cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine proc_addmem(id, size)
      implicit     none
      integer      id, size
      include     'mpif.h'
      include     'mpi_stuff.inc'

      logical      ltrace
      parameter  ( ltrace  = .false. )
      logical      ltrace2
      parameter  ( ltrace2 = .false. )

      if(ltrace2.and.myid.eq.master .or. ltrace2)
     *    write(*,*)'proc_addmem: id, size ', id, size

      proc_meml(id) = proc_meml(id) + size

      if(ltrace.and.myid.eq.master .or. ltrace2) then
          write(*,*) 'proc_addmem: Added workload'
          call proc_dump_wkld()
      end if

      end          ! END: proc_addmem

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  proc_submem:                                                              cc
cc                 Subtracts a certain amount of memory from the storage      cc
cc                 table entry for the given processor. this routines         cc
cc                 pairs with proc_addmem().                                  cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine proc_submem(id, size)
      implicit     none
      integer      id, size
      include     'mpif.h'
      include     'mpi_stuff.inc'

      logical      ltrace
      parameter  ( ltrace  = .false. )
      logical      ltrace2
      parameter  ( ltrace2 = .false. )

      if(ltrace2.and.myid.eq.master .or. ltrace2)
     *    write(*,*)'proc_submem: id, size ', id, size

      proc_meml(id) = proc_meml(id) - size

      if(ltrace.and.myid.eq.master .or. ltrace2) then
          write(*,*) 'proc_addmem: Subtracted memload'
          call proc_dump_wkld()
      end if

      end          ! END: proc_submem

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  proc_dump_wkld:                                                           cc
cc                  Write out the processor workload table and                cc
cc                  minimum loaded processor.                                 cc
cc                  Updated: also updates percentage memory used.             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine proc_dump_wkld()
       implicit   none
       include   'mpif.h'
       include   'mpi_stuff.inc'
       integer    i, totalarena
       !
       integer      mem_return_arena
       external     mem_return_arena

       totalarena = mem_return_arena()

       write(*,99) myid, '     --- Processor workload: ---'
       do i = 0, numprocs-1
          write(*,97) myid, i, proc_wkld(i),100.d0*proc_meml(i)
     *                       /totalarena,
     *                       totalarena-NINT(proc_meml(i))
       end do
       write(*,99) myid, '     Proc w/ minm: ',proc_minwld_id
       write(*,98) myid, '     Minimum load: ',proc_minwld
       write(*,99) myid, '     -- end workload --'

 97      format('[',I4,'] ','     Proc: ',I4,'  Workld: ',F10.2,
c97      format('[',I3,'] ','     Processor: ',I3,'  Workload: ',F16.4,
c    *                      '     %Mem Used: ',F6.1)
     *                      '     %Mem Used: ',F5.1,' Avail:',I11)
 98      format('[',I4,'] ',A,3F16.8)
 99      format('[',I4,'] ',A,3I5)

       return
       end        ! END: proc_dump_wkld

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  proc_init:      Init MPI, and associated variables.                       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine proc_init()
       implicit   none
       include   'mpif.h'
       include   'mpi_stuff.inc'
       integer    i
       integer    OMP_GET_NUM_THREADS
       external   OMP_GET_NUM_THREADS
      logical     ltrace
      parameter ( ltrace  = .false. )

      !
      ! Initialize MPI:
      !
      if(ltrace)write(*,*)'proc_init: Calling MPI_INIT'
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid,     ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      if (numprocs .ge. maxprocs) then
         write(*,*) 'proc_init: numprocs = ', numprocs
         write(*,*) 'proc_init: exceeds maximum allowed: ',maxprocs
         write(*,*) 'proc_init: Increase maxprocs in mpi_stuff.inc '
         write(*,*) 'proc_init: and recompile if possible.'
         write(*,*) 'proc_init: Quitting...'
c        call my_exit('Too many processors')
         stop
      end if
      !  Initialize table of workloads for each processors:
      do i = 0, maxprocs-1
         proc_wkld(i) = 0.d0
      end do
      !  Initialize minimum workload and id w/ least load:
      proc_minwld    = 0.d0
      proc_minwld_id = 0

       numthreads = 1
!$OMP PARALLEL
!$OMP MASTER
!$     numthreads = OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL

       if (myid .eq. master) then
          write(*,*) '  MPI Initialized: numprocs: ',numprocs
          write(*,*) '  Number of threads: ',numthreads
       end if

       return
       end        ! END: proc_init

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  proc_return_wkld:      Return the workload for the calling process.       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       integer function proc_return_wkld()
       implicit   none
       include   'mpif.h'
       include   'mpi_stuff.inc'

       proc_return_wkld = proc_wkld(myid)

       return
       end        ! END: proc_return_wkld

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  proc_return_myid:      Return the id for the calling process.             cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       integer function proc_return_myid()
       implicit   none
       include   'mpif.h'
       include   'mpi_stuff.inc'

       proc_return_myid = myid

       return
       end        ! END: proc_return_myid

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  proc_return_numprocs: Return the number of processes.                     cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       integer function proc_return_numprocs()
       implicit   none
       include   'mpif.h'
       include   'mpi_stuff.inc'

       proc_return_numprocs = numprocs

       return
       end        ! END: proc_return_numprocs

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  listen_for_actions:                                                       cc
cc                      Main driver for request-driven mode.                  cc
cc                      Returns when "DONE" action is received.               cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine listen_for_actions()
      implicit     none
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'action.inc'
      include     'grid_methods.inc'
      include     'glob.inc'
      include     'grid.inc'   ! need for sending gr_error
      include     'param.inc'
      include     'tracers.inc'
      integer  i,              ! generic loop counter
     *         g0,             ! coarse grid identifier
     *         parent,         ! parent of coarse grid (-1 for nil)
     *         action,
     *         gi, gi_tmp,
     *         own, owner,
     *         level, from,
     *         gnx,gny,gnz,
     *         ax,ay,az,
     *         nextproc,
     *         proc_return_wkld,
     *         refinedgrid, deadgrid, deadlevel, numdeadkids, sibling
      integer  j, flev, tmpp

      real(kind=8) tmp_maxchi,   hl,
     *             minx,  miny,  minz,  maxx, maxy, maxz,
     *             minxc, minyc, minzc

      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )
      logical     ltrace3
      parameter ( ltrace3 = .false. )

 99      format('[',I4,'] ',A,4I5)
 98      format('[',I4,'] ',A,I5,A,I5)
 97      format('[',I4,'] ',A,4F15.5)
 96      format('[',I4,'] ',A,F15.5,I5)

         gi = 0
         !
         !  Loop while receiving and executing actions:
         !
 10      continue
         if (ltrace.or.ltrace2.or.ltrace3)
     *                write(*,99) myid,'listen: Waiting for action...'
         call recv_action(action,gi,own,level,parent,
     *                            gnx,gny,gnz,ax,ay,az,from)
         if (ltrace2) write(*,99) myid,'listen: Received:',action, gi
         if (ltrace2) then
            write(*,99) myid,'listen: Received action:',action
            write(*,99) myid,'listen: for grid:       ',gi
            write(*,99) myid,'listen: from processor: ',from
         end if
         flev = level_return_finest(Levelp, maxlev)
         !..............................................................
         if (action .eq. NEWGRID) then
         !..............................................................
            if (ltrace) write(*,99)myid,'listen: NEWGRID ',gi,from
            call send_actionAll(action,gi,own,level,parent,
     *                      gnx,gny,gnz,ax,ay,az)
            !
            ! Determine coordinate bounds:
            !
            if (level.eq.0) then
               hl     = h * refine_factor
               minx   = minx0
               miny   = miny0
               minz   = minz0
               maxx   = maxx0
               maxy   = maxy0
               maxz   = maxz0
            else
               call level_find_bounds(level-1,minx,maxx,miny,
     *                                  maxy,minz,maxz)
               hl     = grid_return_h(level_return_start(level-1))
            end if
            minxc  = minx + (ax-1) * hl
            minyc  = miny + (ay-1) * hl
            minzc  = minz + (az-1) * hl
            if (ltrace2) then
               write(*,96) myid, 'listen:  hl,level: ',hl,level
               write(*,99) myid, 'listen:   ax/y/zc: ',ax,ay,az
               write(*,97) myid, 'listen: minx/y/z : ',minx, miny, minz
               write(*,97) myid, 'listen: minx/y/zc: ',minxc,minyc,minzc
               write(*,97) myid, 'listen: minx/y/z0: ',minx0,miny0,minz0
            end if
            !
            !
            if (myid .ne. own) then
              !
              ! Add new grid to hierarchy and store info:
              !
              if (ltrace) write(*,99)myid,'listen: Just adding:',gi,own
              call grid_newa(gi,own,level,gnx,gny,gnz,minxc,minyc,minzc)
            else
               if (ltrace) write(*,99)myid,'listen: Creating:   ',gi
               !
               ! Create this grid:
               !
               if (ltrace) then
                 write(*,99)myid,'listen:HRun level_test on'
                 call level_test_local()
               end if
               call grid_createA(gi,myid,level,gnx,gny,gnz,
     *                                       minxc,minyc,minzc)
               if (ltrace) write(*,99)myid,'listen: Sending master'
               call send_action(master,action,gi,own,level,parent,
     *                      gnx,gny,gnz,ax,ay,az)
            end if
            if (ltrace2) call level_tree_dump()
            if (ltrace) write(*,99)myid,'listen: Done: NEWGRID'
         !..............................................................
         else if (action .eq. GRIDINIT) then
         !..............................................................
            if (ltrace)write(*,99)myid, 'listen: GRIDINIT from ',from,gi
            call send_actionAll(action,gi,own,level,parent,
     *                      gnx,gny,gnz,ax,ay,az)
            !
            if(ltrace)write(*,99)myid,'listen:enter grid_initialize',gi
            call grid_initialize(gi)
            if (ltrace) then
               write(*,99)myid,'listen:YRun level_test on'
               call level_test_local()
            end if
            if (ltrace)write(*,99)myid, 'listen: GRiDINIT done gi: ',gi
         !..............................................................
         else if (action .eq. MY_BARRIER) then
         !..............................................................
            if (ltrace) write(*,99) myid, 'listen: MY_BARRIER from',from
            call send_actionAll(action,gi,own,level,parent,
     *                      gnx,gny,gnz,ax,ay,az)
            !
            ! Send back acknowledgement
            !
            call MPI_Send(from,1, MPI_INTEGER,
     *                 master, TAG_BARRIER, MPI_COMM_WORLD, ierr)
            if (ltrace) write(*,99) myid, 'listen: MY_BARRIER Finished'
         !..............................................................
         else if (action .eq. DONE) then
         !..............................................................
            if (ltrace) write(*,99) myid, 'listen: Recd DONE from',from
            return
         !..............................................................
         else if (action .eq. DUMP) then
         !..............................................................
            if (ltrace) write(*,99) myid, 'listen: DUMP from',from
            !
            ! Pass info along to next processor:
            !
            !
            ! Dumping current status of processor:
            !
            call level_tree_dump()
            call proc_dump_wkld()
            call mem_stat()
            call level_dump_count()
         !..............................................................
         else if (action .eq. STEP) then
         !..............................................................
            if (ltrace) write(*,99) myid, 'listen: STEP:',level
            call send_actionAll(STEP,gi,0,level,0,0,0,0,0,0,0)
            do i = 1, num_evol_iters
               if (ltrace2) then
                 write(*,99)myid,'listen:LRun level_test on'
                 call level_test_local()
               end if
               if(ltrace)write(*,99) myid, 'listen: STEP: i',i
               call level_iter_local(level)
               if (ltrace2) then
                 write(*,99)myid,'listen:MRun level_test on'
                 call level_test_local()
               end if
               if(ltrace)write(*,99) myid, 'listen: STEP: sync'
               call level_syncbnd_local(level)
               if (ltrace2) then
                 write(*,99)myid,'listen:NRun level_test on'
                 call level_test_local()
               end if
            end do
            if (level .eq. level_return_finest(Levelp,maxlev) .and.
     *          use_mask .gt. 0 ) then
               call level_bhmask_local(level)
               call level_syncbnd_local(level)
            end if
            if (ltrace2) then
                 write(*,99)myid,'listen:ORun level_test on'
                 call level_test_local()
            end if
            if(ltrace)write(*,99) myid, 'listen: STEP: step2'
            call level_step2_local(level)
            if (ltrace) then
                 write(*,99)myid,'listen:PRun level_test on'
                 call level_test_local()
            end if
         !..............................................................
         else if (action .eq. SYNCBND) then
         !..............................................................
            if (ltrace) write(*,99) myid, 'listen: SYNCBND:',level
            call send_actionAll(SYNCBND,gi,0,level,0,0,0,0,0,0,0)
            call level_syncbnd_local(level)
         !..............................................................
         else if (action .eq. INJECT) then
         !..............................................................
            if (ltrace) write(*,99) myid, 'listen: INJECT:',level
            call send_actionAll(INJECT,gi,0,level,0,0,0,0,0,0,0)
            if (ltrace) write(*,99) myid, 'listen: level_inject_local'
            if (ltrace2) call level_tree_dump()
            call level_inject_local(level)
            call level_reset_mask(level-1)
         !..............................................................
         !else if (action .eq. RESETMSK) then
         !..............................................................
            !if (ltrace) write(*,99) myid, 'listen: RESETMSK: ',level
            !call send_actionAll(RESETMSK,gi,0,level,0,0,0,0,0,0,0)
            !call level_reset_mask(level)
         !..............................................................
         else if (action .eq. BOUNDS) then
         !..............................................................
            if (ltrace) write(*,99) myid, 'listen: BOUNDS: ',level
            call send_actionAll(BOUNDS,gi,0,level,0,0,0,0,0,0,0)
            call level_bounds_local(level)
            if (shadow.gt.0) call level_reset_mask(level)
            ! Also need this:
            !call level_syncbnd_local(level)
         !..............................................................
         else if (action .eq. ELLSOLVE) then
         !..............................................................
            if (ltrace) write(*,99) myid, 'listen: ELLSOLVE: ',level
            call send_actionAll(ELLSOLVE,gi,0,level,0,0,0,0,0,0,0)
            if (ltrace) write(*,99) myid,'listen: Calling ell_solve:'
            call ell_solve(level)
            !do j = level_return_finest( Levelp, maxlev ), level+1, -1
               !if(ltrace)write(*,99)myid,'listen: Injecting level: ',j
               !call level_inject_local(j)
               !call level_reset_mask(j-1)
               !call level_get_resid(j-1)  ! just for debugging
            !end do
            if (ltrace) write(*,99) myid,'listen: Done w/ ELLSOLVE'
         !..............................................................
         else if (action .eq. OUTPUT) then
         !..............................................................
            if (ltrace) write(*,99) myid, 'listen: OUTPUT:',level
            call send_actionAll(OUTPUT,gi,0,level,0,0,0,0,0,0,0)
c           if ( (level .eq. 0 .and. shadow .eq. 0) .or.
c    *           (level .eq. 1 .and. shadow .eq. 1) ) then
         if (
     *      (level.eq.0.and.shadow.eq.0 .and. num_anafields.gt.0) .or.
     *      (level.eq.0.and.shadow.eq.1 .and. num_anafields.gt.0
     *                                       .and. allowedl.eq.0) .or.
     *      (level.eq.1.and.shadow.eq.1 .and. num_anafields.gt.0) ) then
               ! If on effective coarse level
               if (ltrace)write(*,*)'listen: Doing analysis on: ',level
               call level_analysis_local(level)
            end if
            ! Surface integrals:
            if(ltrace)write(*,*)'listen:',asf_level,asf_period,level
            if ( level .eq. flev .and. asf_period.gt.0) then
               ! Only compute surfaces after finest level is ready for output
               tmpp   = asf_period
               ! If output requested on finer level than exists,
               ! simply output *every* step of finest level:
               if (asf_level .gt. flev) tmpp = 1
               if (asf_level .lt. flev) then
                  tmpp   = asf_period*2**(flev-asf_level)
               end if
               if(ltrace)write(*,*)'listen:tmp:',flev,tmpp
               if(ltrace)write(*,*)'listen:count',lev_count(flev)
               if (mod( lev_count(flev),tmpp).eq.0) then
                  if(ltrace)write(*,*)'listen:Computing Surface'
                  !call send_actionAll(SURFCOMP,0,myid,0,0,0,0,0,0,0,0)
                  call surface_compute()
                  if(ltrace)write(*,*)'listen:Outputing Surface'
                  ! Only the master outputs these functions
                  !call surface_out(grid_return_time(gi))
               end if
            end if
            !
            !
            ! Tracer particles:  
            if (tracers_period.gt.0) then
               tmpp=tracers_period
               if ( level .eq. 0 .and. 
     .             (mod(lev_count(level)+1,tmpp).eq.0)) then
                  if(ltrace)write(*,*)'listen:: Computing tracers'
                  if(ltrace)write(*,*)'listen::Computing tracers'
                  call tracers_compute(grid_return_time(gi))
                  ! Only master process does this (from level_lib.F):
                  !if(ltrace)write(*,*)'level_apply:Outputing tracers'
                  !call tracers_out(grid_return_time(gi))
                  if(ltrace)write(*,*)'listen:: Done w/ tracers'
               end if
            end if
            !
            call level_output_local(level)
         !..............................................................
         else if (action .eq. COMPUT) then
         !..............................................................
            if (ltrace) write(*,99) myid, 'listen: COMPUT:',level
            call send_actionAll(COMPUT,gi,0,level,0,0,0,0,0,0,0)
            call level_comput_local(level)
            !call level_syncbnd_local(level)
         !..............................................................
         else if (action .eq. LEVINC) then
         !..............................................................
            if (ltrace) write(*,98) myid, 'listen: LEVINC: ',level
            call send_actionAll(LEVINC,gi,0,level,0,0,0,0,0,0,0)
            call level_inc_count(level)
         !..............................................................
         else if (action .eq. SETMSK) then
         !..............................................................
            if (ltrace) write(*,99) myid, 'listen: SETMSK: ',level
            call send_actionAll(SETMSK,gi,0,level,0,0,0,0,0,0,0)
            call level_set_mask(level)
         !..............................................................
         else if (action .eq. BCASTPAR) then
         !..............................................................
            if (ltrace) write(*,99) myid, 'listen: BCASTPAR', from
            call send_actionAll(action,gi,own,level,parent,
     *                      gnx,gny,gnz,ax,ay,az)
            !
            call params_bcast()
            call surface_init()
         !..............................................................
         else if (action .eq. GETFRALL) then
         !..............................................................
            if (ltrace) write(*,99) myid, 'listen: GETFRALL'
            call send_actionAll(action,gi,own,level,parent,
     *                      gnx,gny,gnz,ax,ay,az)
            if (proc_return_wkld().eq.0) then
               maxchi = 0.d0
               if (ltrace2)
     *         write(*,99) myid, 'listen: resetting max for 0 wkld'
            end if
            call MPI_Reduce(maxchi,tmp_maxchi,1,MPI_DOUBLE_PRECISION,
     *                MPI_MAX,0,MPI_COMM_WORLD,ierr)
            !write(*,97) myid, 'listen: maxchi: ',maxchi
            if (ltrace) write(*,99) myid, 'listen: Done MPI_Reduce'
         !..............................................................
         else if (action .eq. GETERROR) then
         !..............................................................
            if (ltrace) write(*,99) myid, 'listen: GETERROR',level
            call send_actionAll(action,gi,own,level,parent,
     *                      gnx,gny,gnz,ax,ay,az)
            if (level .eq. -1)  then
               gi = -1
            else
               gi = level_return_start(level)
            end if
  23        if (grid_return_existence(gi)) then
               if (grid_is_local(gi)) then
                  if (ltrace) write(*,99) myid, 'listen: Error:   ',gi
                  call grid_mask_error(gi) 
                  call load_pointers(gi) 
                  call grid_get_dims(gi,gnx,gny,gnz)
                  call MPI_Send(q(gr_error),gnx*gny*gnz,
     *                 MPI_DOUBLE_PRECISION,
     *                 master, TAG_ERROR_DATA, MPI_COMM_WORLD, ierr)
               end if
               gi = grid_return_sibling(gi)
               goto 23
            end if
            if (ltrace) write(*,99) myid, 'listen: GETERROR done.'
         !..............................................................
         else if (action .eq. FREEGRID) then
         !..............................................................
            if (ltrace) write(*,99) myid, 'listen: Freeing grid ',gi
            call send_actionAll(FREEGRID,gi,0,level,0,0,0,0,0,0,0)
            !
            call level_remove_grid(gi)
            if (ltrace2) call level_tree_dump()
         !..............................................................
         else if (action .eq. BHMASK) then
         !..............................................................
            if (ltrace) write(*,99) myid,'listen: BHMASK: ',level
            call send_actionAll(BHMASK,gi,0,level,0,0,0,0,0,0,0)
            call level_bhmask_local(level)
         !..............................................................
         else if (action .eq. INIT) then
         !..............................................................
            if (ltrace) write(*,99) myid,'listen: INIT: ',level
            call send_actionAll(INIT,gi,0,level,0,0,0,0,0,0,0)
            call level_init_local(level)
         !..............................................................
         else if (action .eq. RESETTIME) then
         !..............................................................
            if (ltrace) write(*,99) myid,'listen: RESETTIME: ',level
            call send_actionAll(RESETTIME,gi,0,level,0,0,0,0,0,0,0)
            call level_resettime_local(level)
         !..............................................................
         else if (action .eq. SAVESTATE) then
         !..............................................................
            if (ltrace) write(*,99) myid,'listen: SAVESTATE: ',level
            call send_actionAll(SAVESTATE,gi,0,level,0,0,0,0,0,0,0)
            call writestate()
         !..............................................................
         else if (action .eq. TRACERINIT) then
         !..............................................................
            if (ltrace) write(*,99) myid,'listen: TRACERINIT: ',level
            call send_actionAll(TRACERINIT,gi,0,level,0,0,0,0,0,0,0)
            call tracers_init()
         !..............................................................
         else if (action .eq. READ_STATE) then
         !..............................................................
            if (ltrace) write(*,99) myid,'listen: READ_STATE: ',level
            call send_actionAll(READ_STATE,gi,0,level,0,0,0,0,0,0,0)
            call readstate()
            write(*,*)'listen: Calling level_modifydata_local()'
            call level_modifydata_local()
            write(*,*)'listen: Done calling level_modifydata_local()'
         !..............................................................
  	      ! Variable timestep modification by Dominic Marcello 2009-08-03
         else if (action .eq. CFL) then
         !..............................................................
            if (ltrace) write(*,99) myid,'listen: CFL: ',level
            call send_actionAll(CFL,gi,0,level,0,0,0,0,0,0,0)
				call comput_cfl()
         !..............................................................
 		   ! --------------end modification by Dominic Marcello 2009-08-03
         else if (action .eq. SURFCOMP) then
         !..............................................................
            if (ltrace) write(*,99) myid,'listen: SURFCOMP: ',level
            call send_actionAll(SURFCOMP,gi,0,level,0,0,0,0,0,0,0)
            call surface_compute()
         !..............................................................
         else if (action .eq. QUIT) then
         !..............................................................
            if (ltrace3) write(*,99) myid,'listen: QUIT: ', from
            !write(*,99) myid,'listen: QUIT: ', from
            call send_actionAll(QUIT,gi,0,level,0,0,0,0,0,0,0)
            call my_exit('Quitting from listen_for_actions.')
         !..............................................................
         else 
         !..............................................................
            write(*,*) 'listen: Received action ',action,' unknown.'
            write(*,*) '                 action   = ',action
            write(*,*) '                 gi       = ',gi
            write(*,*) '                 own      = ',own
            write(*,*) '                 level    = ',level
            write(*,*) '                 parent   = ',parent
            write(*,*) '                 gnx      = ',gnx
            write(*,*) '                 gny      = ',gny
            write(*,*) '                 gnz      = ',gnz
            write(*,*) '                 ax       = ',ax
            write(*,*) '                 ay       = ',ay
            write(*,*) '                 az       = ',az
            write(*,*) 'listen: Quitting.'
            call my_exit('Unknown action.')
         end if
         goto 10
      

      return
      end         ! END: subroutine listen_for_actions

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  proc_getfromall: Called only by master process to                         cc
cc                   get global information.                                  cc
cc                   For this code, it gets the global maxchi                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine  proc_getfromall()
      implicit     none
      real(kind=8)     tmp_maxchi
      include     'action.inc'
      include     'glob.inc'
      include     'mpif.h'
      include     'mpi_stuff.inc'

      logical      ltrace
      parameter  ( ltrace = .false. )

      if (numprocs .eq. 1) return
    
      call send_actionAll(GETFRALL,0,myid,0,0,0,0,0,0,0,0)

      if (ltrace) write(*,*) 'proc_getfromall: MPI_Reduce...',maxchi
      tmp_maxchi = maxchi
      !  The 1 represents the "count"
      !  The 0 represents id of the root, in this case the master
      call MPI_Reduce(tmp_maxchi,maxchi,1,MPI_DOUBLE_PRECISION,MPI_MAX,
     *                0,MPI_COMM_WORLD,ierr)
      if (ltrace) write(*,*) 'proc_getfromall: MPI_Reduce done:',maxchi

      return
      end          ! END: proc_getfromall

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  proc_reconstruct_loads:                                                   cc
cc                 If the processor loads were not stored in the checkpoint,  cc
cc                 then try to reconstruct them by looking at current grid    cc
cc                 structure.                                                 cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine proc_reconstruct_loads()
      implicit     none
      include     'mpif.h'
      include     'mpi_stuff.inc'
      include     'glob.inc'
      include     'largesmall.inc'

      integer      lev, gi, i
      integer      nx,ny,nz
      integer      grid_return_owner, grid_return_sibling
      external     grid_return_owner, grid_return_sibling
      integer      level_return_start
      external     level_return_start
      logical      grid_return_existence
      external     grid_return_existence
      real(kind=8) grid_return_work
      external     grid_return_work
      real(kind=8) loc_proc_wkld(0:numprocs-1),
     *             loc_proc_meml(0:numprocs-1)
      integer      mem_return_meml
      external     mem_return_meml
      logical      doviabroadcast
      parameter  ( doviabroadcast  = .true. )
      logical      ltrace
      parameter  ( ltrace  = .false. )

 98   format('[',I4,']',A,3G12.5)
 99   format('[',I4,']',A,I12)
      if(ltrace.and.myid.eq.master)
     *    write(*,99)myid, 'proc_reconstruct_loads: Enter'

      if (.not.doviabroadcast) then
        write(*,99)myid,'proc_reconstruct_loads: doing locally'
         ! Loop over each level...
         lev = 0
         !    ...on each level, start w/ first grid...
  10     gi  = level_return_start(lev)
         !    ...and update loads for each grid:
  20     if (grid_return_existence(gi)) then
            if(ltrace)write(*,*)myid,'proc_:gi/lev: ',gi,lev
            call grid_get_dims(gi,nx,ny,nz) 
            call proc_addwork( grid_return_owner(gi),
     *                      grid_return_work(gi),
     *                      nx*ny*nz)
            if(ltrace)write(*,*)myid,'proc_:Inc  gi: ',gi
            gi = grid_return_sibling(gi)
            goto 20
         else
            lev = lev + 1
            if(ltrace)write(*,*)myid,'proc_:Inc lev: ',lev
            if (lev .lt. maxlev) goto 10
         end if
         !
      else
         !
         if(ltrace)write(*,99)myid,'proc_reconstruct_loads:doing global'
         if(ltrace)write(*,99)myid,'proc_reconstruct_loads:memload:',
     .                                mem_return_meml()
         if(ltrace)write(*,98)myid,'proc_reconstruct_loads:wkload: ',
     .                                proc_wkld(myid)
         !
         !
         do i = 0, numprocs-1
            if (myid.eq.i) then
               loc_proc_meml(i) = 1.d0*mem_return_meml()
               loc_proc_wkld(i) = proc_wkld(i)
            else
               loc_proc_meml(i) = LARGENUMBER
               loc_proc_wkld(i) = LARGENUMBER
            end if
         end do
         if(ltrace)write(*,98)myid,'proc_reconstruct_loads:locmemload:',
     .                                loc_proc_meml(myid) 
         if(ltrace)write(*,98)myid,'proc_reconstruct_loads:locwkload: ',
     .                                loc_proc_wkld(myid)
         !
         if(ltrace)write(*,99)myid,'proc_reconstruct_loads: proc_meml'
!        call MPI_Reduce(loc_proc_meml(0),proc_meml(0),numprocs,
!    *          MPI_DOUBLE_PRECISION,MPI_MIN,master,MPI_COMM_WORLD,ierr)
         call MPI_Allreduce(loc_proc_meml(0),proc_meml(0),numprocs,
     *          MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
         if(ltrace)write(*,99)myid,'proc_reconstruct_loads: proc_wkld'
!        call MPI_Reduce(loc_proc_wkld(0),proc_wkld(0),numprocs,
!    *          MPI_DOUBLE_PRECISION,MPI_MIN,master,MPI_COMM_WORLD,ierr)
         call MPI_Allreduce(loc_proc_wkld(0),proc_wkld(0),numprocs,
     *          MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
         if(ltrace)write(*,98)myid,'proc_reconstruct_loads:memload:',
     .                                proc_meml(myid)
         if(ltrace)write(*,98)myid,'proc_reconstruct_loads:wkload: ',
     .                                proc_wkld(myid)
      end if

      if(ltrace) then
      !if(ltrace.and.myid.eq.master) then
          write(*,99)myid, 'proc_reconstruct_loads: Done'
          call proc_dump_wkld()
      end if

      end          ! END: proc_reconstruct_lods

