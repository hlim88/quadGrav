cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  mem_return_meml:                                                          cc
cc           Returns the memory load for this processor                       cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function mem_return_meml( )
      implicit    none
      include    'mem.inc'

      mem_return_meml = qptr

      return
      end     ! END: mem_return_meml

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  mem_return_avail:                                                         cc
cc           Returns the size of the available memory in the arena            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function mem_return_avail( )
      implicit    none
      include    'mem.inc'

      mem_return_avail = arena - qptr

      return
      end     ! END: mem_return_avail

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  mem_return_arena:                                                         cc
cc           Returns the size of the arena.                                   cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function mem_return_arena( )
      implicit    none
      include    'mem.inc'

      mem_return_arena = arena

      return
      end     ! END: mem_return_arena

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  mem_stat_brief:                                                           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mem_stat_brief
      implicit none
      include 'grid.inc'
      real(kind=8)   percentused

      percentused = (100.d0/arena)*qptr

      write(*,17)  '                Arena used:     ',qptr,
     *               '/',arena,'  ',percentused,'% (max: ',
     *               (100.d0/arena)*max_qptr,'%)'

 17   format(A,I11,A,I11,A,F4.1,A,F4.1,A)

      return
      end     ! END: subroutine mem_stat_brief

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  mem_stat:                                                                 cc
cc            Print statistics about current memory usage.                    cc
cc            Outputs more information based on parameter "trace_level".      cc
cc                                                                            cc
cc            NB: outputs processor number                                    cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mem_stat
      implicit none
      include 'grid.inc'
      include 'glob.inc'     ! need trace_level variable
      include 'mpif.h'
      include 'mpi_stuff.inc'
      integer  i, sizefree
      real(kind=8)   percentused

      sizefree    = 0.d0

      !
      !  Arena:
      !
      percentused = (100.d0/arena)*qptr
      write(*,9)   myid, 'mem_stat: *****************'
      write(*,17)  myid, 'mem_stat: Arena used:     ',qptr,
     *               '/',arena,'  ',percentused,'% (max: ',
     *               (100.d0/arena)*max_qptr,'%)'
 17   format('[',I4,'] ',A,I11,A,I11,A,F6.2,A,F6.2,A)

      !
      !  Free blocks:
      !
      write(*,16)myid,'mem_stat: # of free blocks:     ',num_freeblocks,
     *              '/',maxfreeblocks,'              (max: ',
     *              max_num_freeblocks,' )'
      do i = 1, num_freeblocks
         percentused = (100.d0/arena)*size_free_blocks(i)
         sizefree    = sizefree + size_free_blocks(i)
         if (trace_level.gt.50) then
         write(*,12) myid,'mem_stat:   Block:       ', i
         write(*,11) myid,'mem_stat:       %used:       ',percentused
         write(*,10) myid,'mem_stat:       FreeBlock:   ',free_blocks(i)
         write(*,10) myid,'mem_stat:       End Block:   ',free_blocks(i)
     *                              + size_free_blocks(i) - 1
         write(*,10) myid,'mem_stat:       SizeBlock:   ',
     *                                size_free_blocks(i)
         end if
      end do
      percentused = (100.d0/arena)*sizefree
      if (num_freeblocks.gt.0) then
         write(*,11)myid,'mem_stat: %Free of blocks:       ',percentused
      end if
      percentused = (100.d0/arena)*size_blockslost
      if (num_blockslost.gt.0) then
         write(*,20)myid,'mem_stat: Blocks lost:          ',
     *         num_blockslost,
     *       '       Size: ',size_blockslost,'  ',percentused,'%'
 20      format('[',I4,'] ',A,I4,A,I9,A,F4.1,A)
      end if

      !
      !  Grid numbers:
      !
      write(*,18)  myid,'mem_stat: # of Grids:           ',numgrids,'/',
     * maxnumgrids,'              (max: ',max_numgrids,' )'
  18  format('[',I4,'] ',A,I5,A,I5,A,I5,A)
      write(*,18)myid, 'mem_stat: # of free grid#s:     ',num_freegrids,
     *    '/',maxnumgrids,'              (max: ',max_num_freegrids,' )'
      if (trace_level.gt.50) then
         do i = 1, num_freegrids
            write(*,12) myid, 'mem_stat:   grid number: ',
     *                  free_grid_numbers(i)
         end do
      end if
      write(*,9)   myid, 'mem_stat: *****************'

   9  format('[',I4,'] ',A)
  10  format('[',I4,'] ',A,I9)
  11  format('[',I4,'] ',A,F6.2)
  12  format('[',I4,'] ',A,I3)
  13  format('[',I4,'] ',A,F6.1,A,F3.0,A)
  14  format('[',I4,'] ',A,I3,A,I3,A)
  15  format('[',I4,'] ',A,I9,A,I9,A)
  16  format('[',I4,'] ',A,I7,A,I7,A,I7,A)

      return
      end     ! END: subroutine mem_stat

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  mem_dealloc:                                                              cc
cc          Frees memory by one of the following:                             cc
cc              --- adding memory to list of free blocks                      cc
cc              --- finding already freed block which is contiguous           cc
cc              --- adding memory back onto arena if at end                   cc
cc              --- disposing of memory if no more storage of free blocks     cc
cc          Doesn't not check that memory wasn't already deallocated.         cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mem_dealloc( pointer, size )
      implicit    none
      integer     pointer, size
      include    'mem.inc'
      integer     endptr, i, smallest, endfreeptr
      integer     myid, proc_return_myid
      external    proc_return_myid

      logical     ltrace
      parameter ( ltrace = .false. )
      
      endptr = pointer + size - 1

 99   format('[',I4,'] ',A,I15)
      if (ltrace) then
         myid = proc_return_myid()
         write(*,99) myid,'mem_dealloc: Freeing memory:'
         write(*,99) myid,'mem_dealloc:    pointer = ', pointer
         write(*,99) myid,'mem_dealloc:    size    = ', size
         write(*,99) myid,'mem_dealloc:    endptr  = ', endptr
         write(*,99) myid,'mem_dealloc:    qptr    = ', qptr
         call mem_stat()
      end if

      !
      ! Check if contiguous with qptr:
      !    if so, then reset qptr
      !    in other words, the block to be freed was last to
      !    be "alloc"ated
      !
      if (qptr .eq. (endptr+1) ) then
         qptr = pointer
         if (ltrace) then
            write(*,99) myid,'mem_dealloc:   Contiguous with qptr'
            write(*,99) myid,'mem_dealloc:   qptr = ',qptr
         end if
         if ( num_freeblocks .gt. 0 ) then
            !
            ! Work backwards looking for free blocks contiq. w/ qptr/arena:
            !
            do i = num_freeblocks, 1, -1
               endfreeptr = free_blocks(i) + size_free_blocks(i) - 1
               If (ltrace) then
                  write(*,99)myid,'mem_dealloc:Checking other blocks.',i
                 write(*,99)myid,'mem_dealloc:endfreeptr+1',endfreeptr+1
               end if
               if ( qptr .eq. (endfreeptr+1) ) then
                  qptr = free_blocks(i)
                  If (ltrace) then
                     write(*,99)myid,'mem_dealloc:Block contig. w/qptr'
                     write(*,99)myid,'mem_dealloc:qptr = ',qptr
                  end if
                  if ( i .eq. num_freeblocks) then
                     ! one less free block (even if num_freeblocks==1)
                     num_freeblocks = num_freeblocks - 1
                  else
                     ! replace block "i" with last block
                     free_blocks(i)      = free_blocks(num_freeblocks)
                     size_free_blocks(i) =
     *                                  size_free_blocks(num_freeblocks)
                     num_freeblocks = num_freeblocks - 1
                  end if
                  If(ltrace)write(*,99)myid,'mem_dealloc:num_freeblock',
     *                                                num_freeblocks
               end if
            end do
         end if
         if (ltrace) then
            write(*,99)myid,'mem_dealloc: Done.'
         end if
         return
      end if
      
      !
      ! If first block to be freed:
      !
      if (num_freeblocks .eq. 0) then
         num_freeblocks      = 1
         max_num_freeblocks  = num_freeblocks
         free_blocks(1)      = pointer
         size_free_blocks(1) = size
         if (ltrace) then
          write(*,99)myid,'mem_dealloc: First block to be freed.'
          write(*,99)myid,'mem_dealloc: num_freeblocks =',num_freeblocks
          write(*,99)myid,'mem_dealloc: free_blocks(1) =',free_blocks(1)
          write(*,99)myid,'mem_dealloc: size_free_blocks(1)=',
     *                              size_free_blocks(1)
          write(*,99)myid,'mem_dealloc: Done.'
         end if
         return
      end if

      !
      ! Check other free blocks to see if contiguous with this block:
      !
      i        = 1
      smallest = 1
  10  endfreeptr = free_blocks(i) + size_free_blocks(i) - 1
      if (free_blocks(i) .eq. (endptr+1)) then
         !
         ! New block contiguous with beginning of an old block:
         !
         free_blocks(i)      = pointer 
         size_free_blocks(i) = size_free_blocks(i) + size
         if (ltrace) then
            write(*,99)myid,'mem_dealloc: Contiguous w/ begin. block',i
            write(*,99)myid,'mem_dealloc:free_blocks(i)=',free_blocks(i)
            write(*,99)myid,'mem_dealloc:size_free_blocks(i)=',
     *                             size_free_blocks(i)
            write(*,99)myid,'mem_dealloc: Done.'
         end if
         return
      else if ( (endfreeptr+1) .eq. pointer) then
         !
         ! New block contiguous with end of an old block:
         !
         size_free_blocks(i) = size_free_blocks(i) + size
         if (ltrace) then
            write(*,99)myid,'mem_dealloc: Contiguous w/ end of block',i
            write(*,99)myid,'mem_dealloc:free_blocks(i)=',free_blocks(i)
            write(*,99)myid,'mem_dealloc:size_free_blocks(i)=',
     *                             size_free_blocks(i)
            write(*,99)myid,'mem_dealloc: Done.'
         end if
         return
      else
         !
         ! Not contiguous with block, look at next one.
         !
         i = i + 1
         if (i.le.num_freeblocks) then
            !
            ! Keep track of smallest sized free block (see below):
            !
            if ( size_free_blocks(smallest) .gt. size_free_blocks(i) ) 
     *           smallest = i
            go to 10
         end if
      end if

      !
      ! No blocks contiguous:
      !
      if (num_freeblocks .lt. maxfreeblocks) then
         if(ltrace)write(*,99)myid,'mem_dealloc:Adding to free list'
         !
         ! Keep track of the free block:
         !
         num_freeblocks                   = num_freeblocks + 1
         if (num_freeblocks.gt.max_num_freeblocks)
     *                 max_num_freeblocks = num_freeblocks
         free_blocks(num_freeblocks)      = pointer
         size_free_blocks(num_freeblocks) = size
         if (ltrace) then
            write(*,99)myid,'mem_dealloc: Adding new free block'
            write(*,99)myid,'mem_dealloc:num_freeblocks=',num_freeblocks
            write(*,99)myid,'mem_dealloc:     free_blocks(n)=',
     *                          free_blocks(num_freeblocks)
            write(*,99)myid,'mem_dealloc:size_free_blocks(n)=',
     *                          size_free_blocks(num_freeblocks)
         end if
      else
         !
         ! Replace the smallest free block
         !      if this block larger than the smallest
         !      else simply lose this memory block
         !
         if (size_free_blocks(smallest) .lt. size) then
            if(ltrace) write(*,99)myid,'mem_dealloc:Replacing ',smallest
            num_blockslost             = num_blockslost + 1
            size_blockslost            =   size_blockslost
     *                                   + size_free_blocks(smallest)
            free_blocks(smallest)      = pointer
            size_free_blocks(smallest) = size
         else
            num_blockslost  = num_blockslost  + 1
            size_blockslost = size_blockslost + size
            if (ltrace) then
               write(*,99)myid,'mem_dealloc: Discarding block, now:'
               write(*,99)myid,'mem_dealloc: num_blockslost:',
     *                                                 num_blockslost
               write(*,99)myid,'mem_dealloc: size_blockslost:',
     *                                                 size_blockslost
            end if
         end if
      end if

      if (ltrace) then
         write(*,99)myid,'mem_dealloc: Done.'
      end if

      return
      end     ! END: mem_dealloc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  mem_alloc:                                                                cc
cc           Try to use a previously freed block.                             cc
cc           If no block big enough, then free in the                         cc
cc           conventional way. If not enough memory, then quit.               cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function mem_alloc( size )
      implicit    none
      integer     size
      include    'mem.inc'
      integer     i,     size_i,
     *            found, size_found, myid
      integer     proc_return_myid
      real(kind=8) percentmemuse 
      external    proc_return_myid
      logical     block_found
      logical     usefreeblocks
      parameter ( usefreeblocks = .true. )
      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )

      block_found = .false.
      myid        = proc_return_myid()

      if (size .lt. 0) then
         write(*,*) 'mem_alloc: Size cannot be negative',size
         call my_exit('mem_alloc: Size cannot be negative')
      end if

      if (ltrace) then
         write(*,*) 'mem_alloc: Allocating memory for block size:',size
         if (.not.usefreeblocks)
     *   write(*,*) 'mem_alloc:      Not using free blocks, only arena'
      end if

      !
      ! If free blocks exist, find smallest block that's big enough
      !
      if ( (num_freeblocks.gt.0) .and. usefreeblocks ) then
         i      = 1
  10     size_i = size_free_blocks(i)
         if (ltrace2)write(*,*) 'mem_alloc:   Looking at free block',i
         if (ltrace2)write(*,*) 'mem_alloc:            with size',size_i
         if (size_i.ge.size) then
            if (ltrace2)write(*,*) 'mem_alloc:    Block large enough...'
            if (block_found) then
               if ( size_i .le. size_found ) then
                  found      = i
                  size_found = size_i
               end if
            else
               block_found = .true.
               found       = i
               size_found  = size_i
            end if
         end if
         i = i + 1
         if (i .le. num_freeblocks) goto 10
      end if

      if (block_found) then
            !
            !  Smallest block which is big enough is grid "found":
            !
            mem_alloc      = free_blocks(found)
            if (ltrace2) then
               write(*,*) 'mem_alloc:       Found free block:'
               write(*,*) 'mem_alloc:     num block found = ',found
               write(*,*) 'mem_alloc: size of block found = ',size_found
               write(*,*) 'mem_alloc:        ptr location = ',mem_alloc
            end if
            if (size_found .eq. size) then
               !
               ! Replace just allocated block with last free block:
               !
               if (ltrace2)write(*,*) 'mem_alloc: free block eliminated'
               free_blocks(found)      = free_blocks(num_freeblocks)
               size_free_blocks(found) =size_free_blocks(num_freeblocks)
               num_freeblocks          = num_freeblocks - 1
            else
               !
               ! Decrease size of the free block:
               !
               free_blocks(found)      = free_blocks(found) + size
               size_free_blocks(found) = size_free_blocks(found) - size
               if (ltrace2) then
                  write(*,*) 'mem_alloc:   size free block decreased:'
                  write(*,*) 'mem_alloc:        free_block(found):',
     *                                               free_blocks(found)
                  write(*,*) 'mem_alloc:   size_free_block(found):',
     *                                          size_free_blocks(found)
               end if
            end if
            if (ltrace2) then
               write(*,*) 'mem_alloc:       Used  free block.'
            end if
      else
         if (ltrace2)write(*,*) 'mem_alloc: Allocating from bulk mem'
         if ((qptr+size).gt.arena) then
              percentmemuse = 100.d0*qptr /arena
              write(*,97)myid,'mem_alloc: Cannot allocate space:'
              write(*,99)myid,'mem_alloc: Size requested: ',size
              write(*,99)myid,'mem_alloc: Arena size:     ',arena
              write(*,99)myid,'mem_alloc:  Memory used:   ',qptr
              write(*,98)myid,'mem_alloc:%Memory used:   ',percentmemuse
              write(*,99)myid,'mem_alloc: myid:           ',myid
              write(*,97)myid,'mem_alloc:If possible,increase arena siz'
              write(*,97)myid,'mem_alloc: in mem.inc, and recompile.'
              call mem_stat
              write(*,99) myid, 'mem_alloc: Quitting.....'
              mem_alloc = -1
              call my_exit('Out of memory.          ')
         else
              mem_alloc = qptr
              qptr      = qptr + size
              if (qptr .gt. max_qptr) max_qptr = qptr
              if (ltrace) write(*,*) 'mem_alloc: qptr = ',qptr
         end if
      end if

      if (ltrace) then
         write(*,*) 'mem_alloc:             mem_alloc = ',mem_alloc
         write(*,*) 'mem_alloc:      mem_alloc + size = ',mem_alloc+size
      end if

      ! If the arena is near the limit for the size of an integer,
      ! it's possible for the arena to overflow and become negative.
      ! This is indicative of running out of memory:
      if (mem_alloc.lt.0) then
         write(*,99) myid, 'mem_alloc: Overflow, out of memory'
         call my_exit('mem_alloc: Out of memory')
      end if

 97   format('[',I4,'] ',A)
 98   format('[',I4,'] ',A,F10.3)
 99   format('[',I4,'] ',A,I15)

      return
      end     ! END: mem_alloc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  mem_writestate:                                                           cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mem_writestate(fname)
      implicit none
      character*(*) fname
      include 'glob.inc'
      include 'grid.inc'
      include 'output.inc'
      include 'mpif.h'
      include 'mpi_stuff.inc'
      include 'mask.inc'
      include 'tracers.inc'
      integer     i, len, ptr, length
      character(5) tmps
      integer     mystringlength
      external    mystringlength
      character(128) name, tmpstr
      !
      !
      !
      logical     ltrace
      parameter ( ltrace = .false. )

      if (ltrace) write(*,*) 'mem_writestate: Begin output to: ',fname

      !
      ! Stuff stored in mem.inc:
      !
      if(ltrace)write(*,*)'mem_: qptr'
      gft_rc=mygft_write_id_int(fname,'qptr',             qptr,    1)
      if(ltrace)write(*,*)'mem_: max_qptr'
      gft_rc=mygft_write_id_int(fname,'max_qptr',         max_qptr,1)
      if(ltrace)write(*,*)'mem_: max_num_freeblocks'
      gft_rc=mygft_write_id_int(fname,'max_num_freeblocks',
     *                                             max_num_freeblocks,1)
      if(ltrace)write(*,*)'mem_: num_freeblocks'
      gft_rc=mygft_write_id_int(fname,'num_freeblocks',   
     *                                             num_freeblocks,1)
      if(ltrace)write(*,*)'mem_: num_blockslost'
      gft_rc=mygft_write_id_int(fname,'num_blockslost',   
     *                                             num_blockslost,1)
      if(ltrace)write(*,*)'mem_: size_blockslost'
      gft_rc=mygft_write_id_int(fname,'size_blockslost',  
     *                                             size_blockslost,1)
      if (num_freeblocks.gt.0) then
      if(ltrace)write(*,*)'mem_: size_free_blocks',num_freeblocks
      gft_rc=mygft_write_id_int(fname,'size_free_blocks',
     *                                             size_free_blocks,
     *                                             num_freeblocks)
      if(ltrace)write(*,*)'mem_: free_blocks',num_freeblocks
      gft_rc=mygft_write_id_int(fname,'free_blocks',
     *                                             free_blocks,
     *                                             num_freeblocks)
      end if

      !****
      ! q
      !****
      if (.not.breakup_q) then
         if(ltrace)write(*,*)'mem_: q ',qptr
         gft_rc=mygft_write_id_float(fname,'q',       q,       qptr    )
      else
         if(ltrace)write(*,*)'mem_: Break q ',qptr
         ptr    = 1
         length = qptr / num_chunks
         do i = 1, num_chunks
            call int2str(i,tmps)
            name = 'q'//tmps
            if (i.eq.num_chunks) length = qptr - (i-1)*length
            if(ltrace)write(*,*)'mem_: q in chunk, i',i,ptr,length,name
            gft_rc=mygft_write_id_float(fname,name,q(ptr),length    )
            ptr = ptr + length
         end do
      end if

      !
      ! Stuff stored in grid.inc:
      !
      if(ltrace)write(*,*)'mem_: grid stuff:'
      gft_rc=mygft_write_id_int(fname,'numgrids',  numgrids,      1)
      gft_rc=mygft_write_id_int(fname,'max_numgrids',     
     *                                             max_numgrids,  1)
      gft_rc=mygft_write_id_int(fname,'num_freegrids',    
     *                                             num_freegrids, 1)
      gft_rc=mygft_write_id_int(fname,'max_num_freegrids',
     *                                             max_num_freegrids,1)
      if (num_freegrids.gt.0)
     *gft_rc=mygft_write_id_int(fname,'free_grid_numbers',
     *                                             free_grid_numbers,
     *                                             num_freegrids)

      gft_rc=mygft_write_id_float(fname,'gr_h',    gr_h,    numgrids)
      gft_rc=mygft_write_id_float(fname,'gr_t',    gr_t,    numgrids)
      gft_rc=mygft_write_id_float(fname,'gr_work', gr_work, numgrids)
      gft_rc=mygft_write_id_float(fname,'gr_maxx', gr_maxx, numgrids)
      gft_rc=mygft_write_id_float(fname,'gr_maxy', gr_maxy, numgrids)
      gft_rc=mygft_write_id_float(fname,'gr_maxz', gr_maxz, numgrids)
      gft_rc=mygft_write_id_float(fname,'gr_minx', gr_minx, numgrids)
      gft_rc=mygft_write_id_float(fname,'gr_miny', gr_miny, numgrids)
      gft_rc=mygft_write_id_float(fname,'gr_minz', gr_minz, numgrids)

      gft_rc=mygft_write_id_int(fname,'gr_own',    gr_own,     numgrids)
      gft_rc=mygft_write_id_int(fname,'gr_count',  gr_count,   numgrids)
      gft_rc=mygft_write_id_int(fname,'gr_level',  gr_level,   numgrids)
      gft_rc=mygft_write_id_int(fname,'gr_parent', gr_parent,  numgrids)
      gft_rc=mygft_write_id_int(fname,'gr_child',  gr_child,   numgrids)
      gft_rc=mygft_write_id_int(fname,'gr_sibling',gr_sibling, numgrids)
      gft_rc=mygft_write_id_int(fname,'gr_size',   gr_size,    numgrids)
      gft_rc=mygft_write_id_int(fname,'gr_nx',     gr_nx,      numgrids)
      gft_rc=mygft_write_id_int(fname,'gr_ny',     gr_ny,      numgrids)
      gft_rc=mygft_write_id_int(fname,'gr_nz',     gr_nz,      numgrids)
      gft_rc=mygft_write_id_int(fname,'gr_x',      gr_x,       numgrids)
      gft_rc=mygft_write_id_int(fname,'gr_y',      gr_y,       numgrids)
      gft_rc=mygft_write_id_int(fname,'gr_z',      gr_z,       numgrids)
      gft_rc=mygft_write_id_int(fname,'gr_alive',  gr_alive,   numgrids)
      gft_rc=mygft_write_id_int(fname,'gr_iter',   gr_iter,    numgrids)

      !
      ! Field pointers from fields.inc:
      !
c     gft_rc=mygft_write_id_int(fname,'gfunc_pointer', gfunc_pointer,
c    *                                    num_gfuncs*maxnumgrids)
      name = 'gfptr'
      len  = mystringlength(name)
      do i = 1, numgrids
         call int2str(i,tmpstr)
         name(len+1:) = tmpstr
         !write(*,*) name
         gft_rc = mygft_write_id_int(fname,name,gfunc_pointer(1,i),
     *                      num_gfuncs)
      end do

      !
      ! Stuff stored in glob.inc:
      !
      gft_rc=mygft_write_id_int(fname,'steps',       steps,         1 )
      gft_rc=mygft_write_id_int(fname,'Levelp',      Levelp(0),  maxlev)
      gft_rc=mygft_write_id_int(fname,'lev_count',  lev_count(0),maxlev)
      gft_rc=mygft_write_id_float(fname,'total_work',total_work,    1 )
      gft_rc=mygft_write_id_float(fname,'maxchi',    maxchi,        1 )
      gft_rc=mygft_write_id_float(fname,'maxchi_t',  maxchi_t,      1 )
      gft_rc=mygft_write_id_int(fname,'maxchi_i',    maxchi_i,      1 )
      gft_rc=mygft_write_id_int(fname,'maxchi_j',    maxchi_j,      1 )
      gft_rc=mygft_write_id_int(fname,'maxchi_k',    maxchi_k,      1 )
      gft_rc=mygft_write_id_int(fname,'maxchi_l',    maxchi_l,      1 )
      gft_rc=mygft_write_id_int(fname,'fileout',     fileout,       1 )

      !
      ! Stuff stored in mpi_stuff.inc:
      !
      gft_rc=mygft_write_id_int(fname,'numprocs',numprocs,       1 )
      if(ltrace)write(*,*)'mem_writestate: proc_wkld'
      gft_rc=mygft_write_id_float(fname,'proc_wkld',
     *                                 proc_wkld(0),numprocs)
      gft_rc=mygft_write_id_float(fname,'proc_meml',
     *                                 proc_meml(0),numprocs)

      !
      ! Stuff stored in mask.inc:
      !
      if(ltrace)write(*,*)'mem_writestate: masking stuff:'
      gft_rc=mygft_write_id_float(fname,'mask_coords',
     *                                 mask_coords,6*max_num_masks)
      gft_rc=mygft_write_id_float(fname,'mask_center',
     *                                 mask_center,3*max_num_masks)
      gft_rc=mygft_write_id_float(fname,'mask_radius',
     *                                 mask_radius,3*max_num_masks)
      gft_rc=mygft_write_id_float(fname,'mask_bbox',
     *                                 mask_bbox,6*max_num_masks)
      gft_rc=mygft_write_id_float(fname,'mask_left',
     *                                 mask_left,3*max_num_masks)
      gft_rc=mygft_write_id_float(fname,'mask_right',
     *                                 mask_right,3*max_num_masks)
      gft_rc=mygft_write_id_float(fname,'local_left',
     *                                 local_left,3*max_num_masks)
      gft_rc=mygft_write_id_float(fname,'local_right',
     *                                 local_right,3*max_num_masks)
      gft_rc=mygft_write_id_float(fname,'transition_time',
     *                                 transition_time,1)

      gft_rc=mygft_write_id_int(fname,'num_masks',     num_masks,     1)
      gft_rc=mygft_write_id_int(fname,'mask_level',    mask_level,    1)
      gft_rc=mygft_write_id_int(fname,'mask_iteration',mask_iteration,1)

      write(*,*) 'mem_writestate: forcedmerge: ',forcedmerge
      write(*,*) 'mem_writestate: bh_true(1): ',bh_true(1)
      write(*,*) 'mem_writestate: bh_true(2): ',bh_true(2)
      gft_rc=mygft_write_id_int(fname,'forcedmerge',   forcedmerge,   1)
      gft_rc=mygft_write_id_int(fname,'bh_true',  bh_true,max_num_masks)

      gft_rc=mygft_write_id_int(fname,'horizon_shape',horizon_shape,2)
      gft_rc=mygft_write_id_int(fname,'horizon_R',  
     *                                 horizon_R,        max_num_masks)
      gft_rc=mygft_write_id_int(fname,'horizon_theta',  
     *                                 horizon_theta,    max_num_masks)
      gft_rc=mygft_write_id_int(fname,'horizon_R_np1',  
     *                                 horizon_R_np1,    max_num_masks)
      gft_rc=mygft_write_id_int(fname,'horizon_theta_np1',  
     *                                 horizon_theta_np1,max_num_masks)
      gft_rc=mygft_write_id_float(fname,'horizon_center',
     *                                 horizon_center,3*max_num_masks)

      !
      ! Stuff stored in tracers.inc:
      !
      if(ltrace)write(*,*)'mem_writestate: tracers stuff:'
      gft_rc=mygft_write_id_float(fname,'tracers_x',
     *                                 tracers_x,tracers_number)
      gft_rc=mygft_write_id_float(fname,'tracers_y',
     *                                 tracers_y,tracers_number)
      gft_rc=mygft_write_id_float(fname,'tracers_z',
     *                                 tracers_z,tracers_number)
      gft_rc=mygft_write_id_float(fname,'tracers_xprev',
     *                                 tracers_xprev,tracers_number)
      gft_rc=mygft_write_id_float(fname,'tracers_yprev',
     *                                 tracers_yprev,tracers_number)
      gft_rc=mygft_write_id_float(fname,'tracers_zprev',
     *                                 tracers_zprev,tracers_number)
      !
      gft_rc=mygft_write_id_float(fname,'tracers_Vx',
     *                                 tracers_Vx,tracers_number)
      gft_rc=mygft_write_id_float(fname,'tracers_Vy',
     *                                 tracers_Vy,tracers_number)
      gft_rc=mygft_write_id_float(fname,'tracers_Vz',
     *                                 tracers_Vz,tracers_number)
      gft_rc=mygft_write_id_float(fname,'tracers_Vxprev',
     *                                 tracers_Vxprev,tracers_number)
      gft_rc=mygft_write_id_float(fname,'tracers_Vyprev',
     *                                 tracers_Vyprev,tracers_number)
      gft_rc=mygft_write_id_float(fname,'tracers_Vzprev',
     *                                 tracers_Vzprev,tracers_number)
      !
      gft_rc=mygft_write_id_float(fname,'tracers_t',
     *                                 tracers_t, 1)

      !
      if (ltrace) write(*,*) 'mem_writestate: End output to: ',fname

      return
      end     ! END: mem_writestate

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                            cc
cc  mem_readstate:                                                            cc
cc                                                                            cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mem_readstate(fname)
      implicit none
      character*(*) fname
      include 'glob.inc'
      include 'grid.inc'
      include 'output.inc'
      include 'mpif.h'
      include 'mpi_stuff.inc'
      include 'largesmall.inc'
      include 'mask.inc'
      include 'tracers.inc'
      integer     i, len, ptr, length, orignumprocs
      integer     mystringlength
      external    mystringlength
      real(kind=8)myl2norm, find_min
      external    myl2norm, find_min
      integer     find_min_index1d
      external    find_min_index1d
      character(5) tmps
      character(128) name, tmpstr
      logical     ltrace
      parameter ( ltrace  = .false. )
      logical     ltrace2
      parameter ( ltrace2 = .false. )

      if(ltrace)write(*,96)myid,'mem_readstate: Begin read from: ',fname

      !
      ! Stuff stored in mem.inc:
      !
      !gft_rc=mygft_read_id_int(fname,'arena',            arena,1)
      if(ltrace2)write(*,99)myid,'mem_readstate: qptr/arena:',qptr,arena
      gft_rc=mygft_read_id_int(fname,'qptr',             qptr,1)
         if (qptr .gt. arena) then
            write(*,99)myid,'mem_readstate: qptr/arena:',qptr,arena
            write(*,97)myid,'mem_readstate: * * * * * * * * * '
            write(*,97)myid,'mem_readstate: Problem: arena is too small'
            write(*,97)myid,'mem_readstate: for this checkpoint.'
            write(*,97)myid,'mem_readstate: Increase and recompile.'
            write(*,97)myid,'mem_readstate: * * * * * * * * * '
            call my_exit('mem_readstate:Arena too smallforcheckpt file')
         end if
      if(ltrace)write(*,99)myid,'mem_readstate: qptr/arena: ',qptr,arena
      gft_rc=mygft_read_id_int(fname,'max_qptr',         max_qptr,1)
      gft_rc=mygft_read_id_int(fname,'max_num_freeblocks',
     *                                           max_num_freeblocks,1)
      gft_rc=mygft_read_id_int(fname,'num_freeblocks',   
     *                                               num_freeblocks,1)
      gft_rc=mygft_read_id_int(fname,'num_blockslost',   
     *                                               num_blockslost,1)
      gft_rc=mygft_read_id_int(fname,'size_blockslost',  
     *                                               size_blockslost,1)

      if (num_freeblocks.gt.0) then
      gft_rc=mygft_read_id_int(fname,'size_free_blocks',
     *                                               size_free_blocks,
     *                                                   num_freeblocks)
      gft_rc=mygft_read_id_int(fname,'free_blocks',free_blocks,
     *                                               num_freeblocks)
      end if

      !****
      ! q
      !****
      if (.not.breakup_q) then
         if(ltrace)write(*,99)myid,'mem_: qptr: ',qptr
         gft_rc=mygft_read_id_float(fname,'q',q, qptr)
      else
         if(ltrace)write(*,99)myid,'mem_: Read in broken q ',qptr
         ptr    = 1
         length = qptr / num_chunks
         do i = 1, num_chunks
            call int2str(i,tmps)
            name = 'q'//tmps
            if (i.eq.num_chunks) length = qptr - (i-1)*length
            if(ltrace)write(*,99)myid,'mem_: q in chunk, i',i,ptr,length
            if(ltrace)write(*,96)myid,'mem_: name:        ',name(1:12)
            gft_rc=mygft_read_id_float(fname,name,q(ptr), length    )
            if(ltrace)write(*,98)myid,'mem_:l2:',myl2norm(q(ptr),length)
            ptr = ptr + length
         end do
      end if
      if(ltrace)write(*,99)myid,'mem_readstate: Done. qptr: ',qptr

      !
      ! Stuff stored in grid.inc:
      !
      gft_rc=mygft_read_id_int(fname,'numgrids',         numgrids,1)
      gft_rc=mygft_read_id_int(fname,'max_numgrids',     max_numgrids,1)
      gft_rc=mygft_read_id_int(fname,'num_freegrids',    
     *                                                num_freegrids,1)
      gft_rc=mygft_read_id_int(fname,'max_num_freegrids',
     *                                             max_num_freegrids,1)
      if (num_freegrids.gt.0)
     *gft_rc=mygft_read_id_int(fname,'free_grid_numbers',
     *                                                free_grid_numbers,
     *                                                    num_freegrids)
      if(ltrace)write(*,99)myid,'mem_readstate: Done grid numbers'

      gft_rc=mygft_read_id_float(fname,'gr_h',   gr_h,    numgrids)
      gft_rc=mygft_read_id_float(fname,'gr_t',   gr_t,    numgrids)
      gft_rc=mygft_read_id_float(fname,'gr_work',gr_work, numgrids)
      gft_rc=mygft_read_id_float(fname,'gr_maxx',gr_maxx, numgrids)
      gft_rc=mygft_read_id_float(fname,'gr_maxy',gr_maxy, numgrids)
      gft_rc=mygft_read_id_float(fname,'gr_maxz',gr_maxz, numgrids)
      gft_rc=mygft_read_id_float(fname,'gr_minx',gr_minx, numgrids)
      gft_rc=mygft_read_id_float(fname,'gr_miny',gr_miny, numgrids)
      gft_rc=mygft_read_id_float(fname,'gr_minz',gr_minz, numgrids)
      if(ltrace)write(*,99)myid,'mem_readstate: Done w/ gr_min/maxes'

      gft_rc=mygft_read_id_int(fname,'gr_own',     gr_own,     numgrids)
      gft_rc=mygft_read_id_int(fname,'gr_count',   gr_count,   numgrids)
      gft_rc=mygft_read_id_int(fname,'gr_level',   gr_level,   numgrids)
      gft_rc=mygft_read_id_int(fname,'gr_parent',  gr_parent,  numgrids)
      gft_rc=mygft_read_id_int(fname,'gr_child',   gr_child,   numgrids)
      gft_rc=mygft_read_id_int(fname,'gr_sibling', gr_sibling, numgrids)
      gft_rc=mygft_read_id_int(fname,'gr_size',    gr_size,    numgrids)
      gft_rc=mygft_read_id_int(fname,'gr_nx',      gr_nx,      numgrids)
      gft_rc=mygft_read_id_int(fname,'gr_ny',      gr_ny,      numgrids)
      gft_rc=mygft_read_id_int(fname,'gr_nz',      gr_nz,      numgrids)
      gft_rc=mygft_read_id_int(fname,'gr_x',       gr_x,       numgrids)
      gft_rc=mygft_read_id_int(fname,'gr_y',       gr_y,       numgrids)
      gft_rc=mygft_read_id_int(fname,'gr_z',       gr_z,       numgrids)
      gft_rc=mygft_read_id_int(fname,'gr_alive',   gr_alive,   numgrids)
      gft_rc=mygft_read_id_int(fname,'gr_iter',    gr_iter,    numgrids)
      if(ltrace)write(*,99)myid,'mem_readstate: Done w/ grids'

      !
      ! Field pointers from fields.inc:
      !
c     gft_rc=mygft_read_id_int(fname,'gfunc_pointer', gfunc_pointer,
c    *                                    num_gfuncs*maxnumgrids)
      name = 'gfptr'
      len  = mystringlength(name)
      do i = 1, numgrids
         call int2str(i,tmpstr)
         name(len+1:) = tmpstr
      gft_rc=mygft_read_id_int(fname,name,gfunc_pointer(1,i),num_gfuncs)
      end do
      if(ltrace)write(*,99)myid,'mem_readstate: Done w/ w/ field ptrs'

      !
      ! Stuff stored in glob.inc:
      !
      gft_rc=mygft_read_id_int(fname,'steps',        steps,         1 )
      gft_rc=mygft_read_id_int(fname,'Levelp',      Levelp(0),   maxlev)
      gft_rc=mygft_read_id_int(fname,'lev_count',   lev_count(0),maxlev)
      gft_rc=mygft_read_id_float(fname,'total_work', total_work,    1 )
      gft_rc=mygft_read_id_float(fname,'maxchi',     maxchi,        1 )
      gft_rc=mygft_read_id_float(fname,'maxchi_t',   maxchi_t,      1 )
      gft_rc=mygft_read_id_int(fname,'maxchi_i',     maxchi_i,      1 )
      gft_rc=mygft_read_id_int(fname,'maxchi_j',     maxchi_j,      1 )
      gft_rc=mygft_read_id_int(fname,'maxchi_k',     maxchi_k,      1 )
      gft_rc=mygft_read_id_int(fname,'maxchi_l',     maxchi_l,      1 )
      gft_rc=mygft_read_id_int(fname,'fileout',      fileout,       1 )
      if(ltrace)write(*,99)myid,'mem_readstate: Done w/ w/ glob.inc'


      if (.true.) then
         if (ltrace) then
            write(*,96)myid,'mem_readstate:           * * * * * '
            write(*,96)myid,'mem_readstate: To allow for an increase to'
            write(*,96)myid,'mem_readstate: the number of processes for'
            write(*,96)myid,'mem_readstate: a restart from checkpoints,'
            write(*,96)myid,'mem_readstate: we need to broadcast some'
            write(*,96)myid,'mem_readstate: info...'
            write(*,96)myid,'mem_readstate: (if the code fails here,'
            write(*,96)myid,'mem_readstate:check for explanation above)'
            write(*,96)myid,'mem_readstate: Calling MPI_BARRIER'
         end if
         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         if(ltrace)write(*,96)myid,'mem_readstate:...back from Barrier'
         if(ltrace)write(*,96)myid,'mem_readstate:...AMR hierarchy...'
         if(ltrace2)write(*,99)myid,'numgrids',numgrids
         call MPI_BCAST(numgrids,     1, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
         if(ltrace)write(*,99)myid,'numgrids',numgrids
         !
         if(ltrace2)write(*,99)myid,'max_numgrids',max_numgrids
         call MPI_BCAST(max_numgrids, 1, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
         if(ltrace)write(*,99)myid,'max_numgrids',max_numgrids
         !
         if(ltrace2)write(*,99)myid,'num_freegrids',num_freegrids
         call MPI_BCAST(num_freegrids,1, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
         if(ltrace)write(*,99)myid,'num_freegrids',num_freegrids
         !
         if(ltrace2)write(*,99)myid,'max_num_freegris',max_num_freegrids
         call MPI_BCAST(max_num_freegrids,1, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
         if(ltrace)write(*,99)myid,'max_num_freegrids',max_num_freegrids
         !
         if (num_freegrids.gt.0)
     *   call MPI_BCAST(free_grid_numbers,num_freegrids, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
         !
         if (numgrids.gt.0) then
            if (ltrace2) call level_tree_dump()
            if(ltrace2)write(*,98)myid,'gr_h():',
     *                            gr_h(1),gr_h(2),gr_h(3)
         call MPI_BCAST(gr_h,            numgrids, MPI_DOUBLE_PRECISION,
     &                                  master, MPI_COMM_WORLD, ierr)
            if(ltrace)write(*,98)myid,'gr_h():',
     *                            gr_h(1),gr_h(2),gr_h(3)
            if(ltrace2)write(*,98)myid,'gr_t():',
     *                            gr_t(1),gr_t(2),gr_t(3)
         call MPI_BCAST(gr_t,            numgrids, MPI_DOUBLE_PRECISION,
     &                                  master, MPI_COMM_WORLD, ierr)
            if(ltrace2)write(*,98)myid,'gr_t():',
     *                            gr_t(1),gr_t(2),gr_t(3)
         call MPI_BCAST(gr_work,         numgrids, MPI_DOUBLE_PRECISION,
     &                                  master, MPI_COMM_WORLD, ierr)
            if(ltrace2)write(*,98)myid,'gr_maxx():',
     *                            gr_maxx(1),gr_maxx(2),gr_maxx(3)
         call MPI_BCAST(gr_maxx,         numgrids, MPI_DOUBLE_PRECISION,
     &                                  master, MPI_COMM_WORLD, ierr)
            if(ltrace2)write(*,98)myid,'gr_maxx():',
     *                            gr_maxx(1),gr_maxx(2),gr_maxx(3)
            if(ltrace2)write(*,98)myid,'gr_maxy():',
     *                            gr_maxy(1),gr_maxy(2),gr_maxy(3)
         call MPI_BCAST(gr_maxy,         numgrids, MPI_DOUBLE_PRECISION,
     &                                  master, MPI_COMM_WORLD, ierr)
            if(ltrace2)write(*,98)myid,'gr_maxy():',
     *                            gr_maxy(1),gr_maxy(2),gr_maxy(3)
            if(ltrace2)write(*,98)myid,'gr_maxz():',
     *                            gr_maxz(1),gr_maxz(2),gr_maxz(3)
         call MPI_BCAST(gr_maxz,         numgrids, MPI_DOUBLE_PRECISION,
     &                                  master, MPI_COMM_WORLD, ierr)
            if(ltrace2)write(*,98)myid,'gr_maxz():',
     *                            gr_maxz(1),gr_maxz(2),gr_maxz(3)
            if(ltrace2)write(*,98)myid,'gr_minx():',
     *                            gr_minx(1),gr_minx(2),gr_minx(3)
         call MPI_BCAST(gr_minx,         numgrids, MPI_DOUBLE_PRECISION,
     &                                  master, MPI_COMM_WORLD, ierr)
            if(ltrace2)write(*,98)myid,'gr_minx():',
     *                            gr_minx(1),gr_minx(2),gr_minx(3)
            if(ltrace2)write(*,98)myid,'gr_miny():',
     *                            gr_miny(1),gr_miny(2),gr_miny(3)
         call MPI_BCAST(gr_miny,         numgrids, MPI_DOUBLE_PRECISION,
     &                                  master, MPI_COMM_WORLD, ierr)
            if(ltrace2)write(*,98)myid,'gr_miny():',
     *                            gr_miny(1),gr_miny(2),gr_miny(3)
            if(ltrace2)write(*,98)myid,'gr_minz():',
     *                            gr_minz(1),gr_minz(2),gr_minz(3)
         call MPI_BCAST(gr_minz,         numgrids, MPI_DOUBLE_PRECISION,
     &                                  master, MPI_COMM_WORLD, ierr)
            if(ltrace2)write(*,98)myid,'gr_minz():',
     *                            gr_minz(1),gr_minz(2),gr_minz(3)
         !
            if(ltrace2)write(*,99)myid,'gr_own():',
     *                            gr_own(1),gr_own(2),gr_own(3)
         call MPI_BCAST(gr_own,           numgrids, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
            if(ltrace2)write(*,99)myid,'gr_own():',
     *                            gr_own(1),gr_own(2),gr_own(3)
         call MPI_BCAST(gr_count,         numgrids, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(gr_level,         numgrids, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(gr_parent,         numgrids, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(gr_child,          numgrids, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(gr_sibling,        numgrids, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(gr_size,           numgrids, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(gr_nx,             numgrids, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(gr_ny,             numgrids, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(gr_nz,             numgrids, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(gr_alive,          numgrids, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(gr_iter,           numgrids, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
            call level_tree_dump()
         end if
         if(ltrace2)write(*,99)myid,'steps            ',steps
         call MPI_BCAST(steps,            1, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
         if(ltrace2)write(*,99)myid,'steps            ',steps
         if(maxlev.gt.0) then
            if(ltrace2)write(*,99)myid,'Levelp():',
     *                            Levelp(0),Levelp(1),Levelp(2)
            call MPI_BCAST(Levelp(0),   maxlev, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
            if(ltrace)write(*,99)myid,'Levelp():',
     *                            Levelp(0),Levelp(1),Levelp(2)
            if(ltrace2)write(*,99)myid,'lev_count():',
     *                            lev_count(0),lev_count(1),lev_count(2)
            call MPI_BCAST(lev_count(0),maxlev, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
            if(ltrace2)write(*,99)myid,'lev_count():',
     *                            lev_count(0),lev_count(1),lev_count(2)
         end if
         if (ltrace) then
            call level_tree_dump()
            write(*,96)myid,'mem_readstate: Done w/ broadcasting.'
         end if
      end if
      !
      ! Stuff stored in mpi_stuff.inc...
      !
      gft_rc=mygft_read_id_int(fname,'numprocs', orignumprocs,1 )
      if (gft_rc.lt.1) orignumprocs=numprocs
      if(ltrace)write(*,97)myid,'mem_readstate: proc_wkld'
      gft_rc=mygft_read_id_float(fname,'proc_wkld',
     *                                 proc_wkld(0),orignumprocs)
!     gft_rc=mygft_read_id_float(fname,'proc_meml',
!    *                                 proc_meml(0),numprocs)
      ! ...if these were not found, then try to reconstruct...
      !if (proc_wkld(0) .le. SMALLNUMBER) then
      if (.true.) then
         if(ltrace)write(*,97)myid,'mem_readstate: reconstruct_loads'
         call proc_reconstruct_loads()
      else if(ltrace) then
         write(*,97)myid,'mem_readstate: Workloads already loaded'
         write(*,98)myid,'mem_readstate: proc_wkld(0)=',proc_wkld(0)
      end if
      !if (ltrace.and.myid.eq.master) call proc_dump_wkld()
      ! ...and find the one w/ least load:
      !     (NB: subtract one because of indexing from 0)
      proc_minwld    = find_min(        proc_wkld(0),numprocs)
      proc_minwld_id = find_min_index1D(proc_wkld(0),numprocs)-1
      if (ltrace) then
         write(*,98)myid,'mem_readstate:proc_minwld=',proc_minwld
         write(*,99)myid,'mem_readstate:proc_minwld_id=',proc_minwld_id
      end if

      !
      ! Stuff stored in mask.inc:
      !
      if(ltrace)write(*,97)myid,'mem_readstate: masking stuff:'
      gft_rc=mygft_read_id_float(fname,'mask_coords',
     *                                 mask_coords,6*max_num_masks)
      call MPI_BCAST(mask_coords, 6*max_num_masks, MPI_DOUBLE_PRECISION,
     &                                  master, MPI_COMM_WORLD, ierr)
      !
      gft_rc=mygft_read_id_float(fname,'mask_center',
     *                                 mask_center,3*max_num_masks)
      call MPI_BCAST(mask_center, 3*max_num_masks, MPI_DOUBLE_PRECISION,
     &                                  master, MPI_COMM_WORLD, ierr)
      !
      gft_rc=mygft_read_id_float(fname,'mask_radius',
     *                                 mask_radius,3*max_num_masks)
      call MPI_BCAST(mask_radius, 3*max_num_masks, MPI_DOUBLE_PRECISION,
     &                                  master, MPI_COMM_WORLD, ierr)
      !
      gft_rc=mygft_read_id_float(fname,'mask_bbox',
     *                                 mask_bbox,6*max_num_masks)
      call MPI_BCAST(mask_bbox,   6*max_num_masks, MPI_DOUBLE_PRECISION,
     &                                  master, MPI_COMM_WORLD, ierr)
      !
      gft_rc=mygft_read_id_float(fname,'mask_left',
     *                                 mask_left,3*max_num_masks)
      call MPI_BCAST(mask_left,   3*max_num_masks, MPI_DOUBLE_PRECISION,
     &                                  master, MPI_COMM_WORLD, ierr)
      !
      gft_rc=mygft_read_id_float(fname,'mask_right',
     *                                 mask_right,3*max_num_masks)
      call MPI_BCAST(mask_right,  3*max_num_masks, MPI_DOUBLE_PRECISION,
     &                                  master, MPI_COMM_WORLD, ierr)
      !
      gft_rc=mygft_read_id_float(fname,'local_left',
     *                                 local_left,3*max_num_masks)
      call MPI_BCAST(local_left,  3*max_num_masks, MPI_DOUBLE_PRECISION,
     &                                  master, MPI_COMM_WORLD, ierr)
      !
      gft_rc=mygft_read_id_float(fname,'local_right',
     *                                 local_right,3*max_num_masks)
      call MPI_BCAST(local_right, 3*max_num_masks, MPI_DOUBLE_PRECISION,
     &                                  master, MPI_COMM_WORLD, ierr)
      !
      gft_rc=mygft_read_id_float(fname,'transition_time',
     *                                 transition_time,1)
      call MPI_BCAST(transition_time,           1, MPI_DOUBLE_PRECISION,
     &                                  master, MPI_COMM_WORLD, ierr)
      !

      gft_rc=mygft_read_id_int(fname,'num_masks',     num_masks,     1)
      call MPI_BCAST(num_masks,                 1, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
      !
      gft_rc=mygft_read_id_int(fname,'mask_level',    mask_level,    1)
      call MPI_BCAST(mask_level,                1, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
      !
      gft_rc=mygft_read_id_int(fname,'mask_iteration',mask_iteration,1)
      call MPI_BCAST(mask_iteration,            1, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
      !

      gft_rc=mygft_read_id_int(fname,'forcedmerge',   forcedmerge,   1)
      call MPI_BCAST(forcedmerge,               1, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
      !
      gft_rc=mygft_read_id_int(fname,'bh_true',  bh_true,max_num_masks)
      call MPI_BCAST(bh_true,                   1, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
      !
      if (ltrace) then
        write(*,99)myid,'mem_readstate: forcedmerge = ',forcedmerge
        write(*,99)myid,'mem_readstate: bh_true(1)  = ',bh_true(1)
        write(*,99)myid,'mem_readstate: bh_true(2)  = ',bh_true(2)
      end if

      gft_rc=mygft_read_id_int(fname,'horizon_shape',horizon_shape,2)
      call MPI_BCAST(horizon_shape,             2, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
      !
      gft_rc=mygft_read_id_int(fname,'horizon_R',  
     *                                 horizon_R,        max_num_masks)
      call MPI_BCAST(horizon_R,     max_num_masks, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
      !
      gft_rc=mygft_read_id_int(fname,'horizon_theta',  
     *                                 horizon_theta,    max_num_masks)
      call MPI_BCAST(horizon_theta, max_num_masks, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
      !
      gft_rc=mygft_read_id_int(fname,'horizon_R_np1',  
     *                                 horizon_R_np1,    max_num_masks)
      call MPI_BCAST(horizon_R_np1, max_num_masks, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
      !
      gft_rc=mygft_read_id_int(fname,'horizon_theta_np1',  
     *                                 horizon_theta_np1,max_num_masks)
      call MPI_BCAST(horizon_theta_np1, max_num_masks, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
      !
      gft_rc=mygft_read_id_float(fname,'horizon_center',
     *                                 horizon_center,3*max_num_masks)
      call MPI_BCAST(horizon_center,3*max_num_masks, MPI_INTEGER,
     &                                  master, MPI_COMM_WORLD, ierr)
      !

      !
      ! Stuff stored in tracers.inc:
      !
!     common  / tracers /  tracers_x,
!    *                     tracers_y,
!    *                     tracers_z,
!    *                     tracers_xprev,
!    *                     tracers_yprev,
!    *                     tracers_zprev,
!    *                     tracers_t
      if(ltrace)write(*,*)'mem_readstate: tracers stuff:'
      gft_rc=mygft_read_id_float(fname,'tracers_x',
     *                                 tracers_x,tracers_number)
      gft_rc=mygft_read_id_float(fname,'tracers_y',
     *                                 tracers_y,tracers_number)
      gft_rc=mygft_read_id_float(fname,'tracers_z',
     *                                 tracers_z,tracers_number)
      gft_rc=mygft_read_id_float(fname,'tracers_xprev',
     *                                 tracers_xprev,tracers_number)
      gft_rc=mygft_read_id_float(fname,'tracers_yprev',
     *                                 tracers_yprev,tracers_number)
      gft_rc=mygft_read_id_float(fname,'tracers_zprev',
     *                                 tracers_zprev,tracers_number)
      gft_rc=mygft_read_id_float(fname,'tracers_t',
     *                                 tracers_t, 1)


      !
      !
      !
      if (ltrace2.and.myid.eq.master) call level_tree_dump()
      if(ltrace)write(*,96)myid,'mem_readstate: End   read from: ',fname

 96   format('[',I4,'] ',A,A)
 97   format('[',I4,'] ',A)
 98   format('[',I4,'] ',A,3F10.3)
 99   format('[',I4,'] ',A,3I15)
      return
      end     ! END: mem_readstate

