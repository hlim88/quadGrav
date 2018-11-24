      ! "Dynamic" memory:
      
      
      
      
      
      
      
      ! size of "dynamic" memory
      integer     arena
      parameter ( arena          = 220000000 )
      
      integer     maxfreeblocks
      parameter ( maxfreeblocks  = 4000        )
      
      real(kind=8)        q(arena)
      common      / mem_common /
     *              q
      
      integer       free_blocks(maxfreeblocks)
      common      / mem_common /
     *              free_blocks
      
      integer       size_free_blocks(maxfreeblocks) 
      common      / mem_common /
     *              size_free_blocks
      
      
      
      
      ! it will write "q" in that number of chunks.
      
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
