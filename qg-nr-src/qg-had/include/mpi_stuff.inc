c----------------------------------------------------------------------
c     Variables required for MPI and such:
c     Must include mpif.h first (      include     'mpif.h'    )
c      else get an error about MPI_STATUS_SIZE
c
c----------------------------------------------------------------------
      integer     num_int_params,       num_dbl_params
      parameter ( num_int_params = 200, num_dbl_params = 300 )

      ! temporary storage for bcast params
      integer     tmp_ints(num_int_params)
      ! temporary storage for bcast params
      real(kind=8)      tmp_dbls(num_dbl_params)

      integer     ierr,          ! return code from MPI calls
     *            myid,          ! non-negative ID for this processor
     *            numprocs,      ! number of processors for this run
     *            proc_minwld_id,! id of least loaded processor
     *            numthreads     ! number of threads for this run

      common   /  mpistuff_ints /
     *            tmp_ints,
     *            ierr,
     *            myid,
     *            numprocs,
     *            proc_minwld_id

      ! maximum number of allowed procs
      integer     maxprocs
      parameter ( maxprocs      =  4096 )

      ! null processor --- no communication
      integer     NULL_PROC
      parameter ( NULL_PROC     =  -1 )

      ! ID for master processor
      integer     master
      parameter ( master        =  0 )

c     ! workload for each processor and minimum current workload
      real(kind=8)      proc_wkld(0:maxprocs-1), proc_minwld

c     ! memory load  for each processor
      real(kind=8)      proc_meml(0:maxprocs-1)

      common   /  mpistuff_dbls /
     *            tmp_dbls,
     *            proc_wkld,
     *            proc_meml,
     *            proc_minwld,
     *            numthreads

         !
         !   Format of a grid message:
         !         gridmessage( 1) = action
         !         gridmessage( 2) = grid number
         !         gridmessage( 3) = owner process
         !         gridmessage( 4) = level of grid
         !         gridmessage( 5) = parent of grid
         !         gridmessage( 6) = nx
         !         gridmessage( 7) = ny
         !         gridmessage( 8) = nz
         !         gridmessage( 9) = ax
         !         gridmessage(10) = ay
         !         gridmessage(11) = az
         !

      ! msg length for grid actions
      integer     mlength
      parameter ( mlength = 11 )
      ! actual msg which gets sent
      integer     gridmsg(mlength),
     *            status(MPI_STATUS_SIZE)

      !
      !  Types of tags for the MPI communications:
      !
      integer     TAG_ACTION,
     *            TAG_BND_DATA,
     *            TAG_INJ_DATA,
     *            TAG_INIT_DATA,
     *            TAG_SYNC_DATA,
     *            TAG_ERROR_DATA,
     *            TAG_CGC_DATA,
     *            TAG_PST_DATA,
     *            TAG_GRID_AWARE,
     *            TAG_BARRIER,
     *            TAG_SURFACE,
     *            TAG_CHKPT,
     *            TAG_GATHER
      parameter ( TAG_ACTION     = 0,
     *            TAG_BND_DATA   = 1,
     *            TAG_INJ_DATA   = 2,
     *            TAG_INIT_DATA  = 3,
     *            TAG_SYNC_DATA  = 4,
     *            TAG_ERROR_DATA = 5,
     *            TAG_CGC_DATA   = 6,
     *            TAG_PST_DATA   = 7,
     *            TAG_GRID_AWARE = 8,
     *            TAG_BARRIER    = 9,
     *            TAG_SURFACE    = 10,
     *            TAG_CHKPT      = 11,
     *            TAG_GATHER     = 12
     *            )
