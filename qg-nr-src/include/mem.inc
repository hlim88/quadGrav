      ! "Dynamic" memory:
      !
      !                            NB: Changes here demand similar
      !                                changes to the mem_writestate()
      !                                           and mem_readstate()
      !                                routines.                       
      !
      !

      ! size of "dynamic" memory
      integer     arena

!     arena          =  10 000 000       --size nlsm:  82   MB
!     arena          =  25 000 000       --size nlsm: 201   MB
!     arena          =  30 000 000       --size nlsm: 241   MB
!     arena          =  50 000 000       --size nlsm: 402   MB
!     arena          =  90 000 000       --size nlsm: 721   MB
!     arena          = 110 000 000       --size nlsm: 881   MB
!     arena          = 170 000 000       --size nlsm  1.36  GB
!     arena          = 200 000 000       --size nlsm  1.6   GB
!     arena          = 220 000 000       --size nlsm  1.76  GB
!     arena          = 820 000 000       --size nlsm  6.58  GB
!     arena          = 2 147 000 000     --size nlsm 17.2   GB
!     arena          = 7 147 000 000     --size nlsm 57.2   GB

!     parameter ( arena          =  10 000 000 )
!     parameter ( arena          =  25 000 000 )
!     parameter ( arena          =  30 000 000 )
!     parameter ( arena          =  50 000 000 )
!      parameter ( arena          =  90 000 000 )
!     parameter ( arena          = 110 000 000 )
!     parameter ( arena          = 170 000 000 )
!     parameter ( arena          = 200 000 000 )
      parameter ( arena          = 220 000 000 )
!     parameter ( arena          = 820 000 000 )
!     2^63 = 9223372036854775808   -integer-size 64
!     2^31 =          2147483648   
!     parameter ( arena          = 2 147 000 000 )
!     parameter ( arena          = 7 147 000 000 )


      ! total # of freed blocks any more are lost
      integer     maxfreeblocks
      parameter ( maxfreeblocks  = 4000        )

      ! actual memory
      real(kind=8)        q(arena)
      common      / mem_common /
     *              q

      ! pointers to free blocks
      integer       free_blocks(maxfreeblocks)
      common      / mem_common /
     *              free_blocks

      ! pointers to free blocks
      integer       size_free_blocks(maxfreeblocks) 
      common      / mem_common /
     *              size_free_blocks

!   qptr,                           ! pointer into dynamic mem
!   max_qptr,                       ! max pointer at one time
!   num_freeblocks,                 ! # of free blocks available
!   max_num_freeblocks,             ! max # blocks at 1 time
!   num_blockslost,                 ! # of free blocks lost
!   size_blockslost                 ! amount of memory lost


      !
      ! User options for how to checkpoint
      ! the large memory arena on the F77 side
      ! If you set breakup_q to true, then
      ! it will write "q" in that number of chunks.
      !
      logical       breakup_q
      parameter (   breakup_q = .true. )
      integer       num_chunks
      parameter (   num_chunks = 100 )

      integer       qptr,
     *              max_qptr,
     *              num_freeblocks,
     *              max_num_freeblocks,
     *              num_blockslost,
     *              size_blockslost
      common      / mem_common /
     *              qptr,                       
     *              max_qptr,                   
     *              num_freeblocks,              
     *              max_num_freeblocks,           
     *              num_blockslost,                
     *              size_blockslost                 
